import argparse
import torch
import platform
import re
import logging
import numpy as np
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional
from Bio import SeqIO
import pandas as pd
from scipy.spatial.distance import cdist
import time
import os
from multiprocessing import Pool

logger = logging.getLogger(__name__)

def get_device() -> torch.device:
    """Determine the appropriate device for computation."""
    if torch.cuda.is_available():
        return torch.device('cuda')
    elif torch.backends.mps.is_available() and platform.system() == 'Darwin':
        return torch.device('mps')
    return torch.device('cpu')

class EmbeddingError(Exception):
    """Raised when there's an error during embedding generation."""
    pass

def normalize_l2(x: np.ndarray) -> np.ndarray:
    """Normalize vector using L2 normalization."""
    norm = np.linalg.norm(x)
    return x / (norm if norm > 0 else 1)

class ProteinEmbedder(ABC):
    """Abstract base class for protein embedders."""

    def __init__(self):
        self.device = get_device()
        logger.info(f"Initializing {self.__class__.__name__} using device: {self.device}")

    @property
    @abstractmethod
    def vector_size(self) -> int:
        """Return the size of the generated embedding vector."""
        pass

    @abstractmethod
    def get_embedding(self, sequence: str) -> List[float]:
        """Generate embedding for a protein sequence."""
        pass

    @abstractmethod
    def get_batch_embeddings(self, sequences: List[str], batch_size: int = 8) -> List[List[float]]:
        """Generate embeddings for multiple sequences in batches."""
        pass

    def validate_embedding(self, embedding: List[float]) -> List[float]:
        """Validate embeddings to ensure they are well-formed."""
        if not embedding:
            raise EmbeddingError("Generated embedding is empty")
        if not all(isinstance(x, (float, np.floating)) for x in embedding):
            # Convert to float if needed
            embedding = [float(x) for x in embedding]
        if any(np.isnan(x) or np.isinf(x) for x in embedding):
            raise EmbeddingError("NaN or Inf values in embedding")
        return embedding

class ESM2Embedder(ProteinEmbedder):
    """Optimized ESM2 embedder with batch processing."""

    def __init__(self, model_name="facebook/esm2_t12_35M_UR50D"):
        self.device = get_device()
        logger.info(f"Initializing ESM2Embedder using device: {self.device}")

        try:
            from transformers import AutoTokenizer, AutoModel
            self.model_name = model_name
            self.tokenizer = AutoTokenizer.from_pretrained(model_name)
            self.model = AutoModel.from_pretrained(model_name)
            self.model = self.model.to(self.device)
            self.model.eval()
            self._vector_size = self.model.config.hidden_size

        except Exception as e:
            logger.error(f"Failed to initialize ESM2 model: {str(e)}")
            raise EmbeddingError(f"Failed to initialize ESM2 model: {str(e)}")

    @property
    def vector_size(self) -> int:
        return self._vector_size

    def get_embedding(self, sequence: str, max_length: int = 1024) -> List[float]:
        """Handle long sequences by chunking and averaging embeddings."""
        seq = re.sub(r'[^A-Z]', '', sequence.upper())
        seq = re.sub(r'[UZOB]', 'X', seq)
        chunks = [seq[i:i + max_length - 2] for i in range(0, len(seq), max_length - 2)]
        chunk_embeddings = []
        for chunk in chunks:
            emb = self._embed_single(chunk, max_length)
            chunk_embeddings.append(emb)
        mean_embedding = np.mean(np.array(chunk_embeddings), axis=0)
        normalized_embedding = normalize_l2(mean_embedding)
        return self.validate_embedding(normalized_embedding.tolist())

    def _embed_single(self, seq: str, max_length: int) -> np.ndarray:
        inputs = self.tokenizer(
            seq,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=max_length,
            add_special_tokens=True
        )
        inputs = {key: val.to(self.device) for key, val in inputs.items()}
        with torch.no_grad():
            outputs = self.model(**inputs)
        last_hidden_state = outputs.last_hidden_state
        attention_mask = inputs['attention_mask']
        valid_positions = attention_mask[0].bool()
        sequence_embeddings = last_hidden_state[0][valid_positions]
        mean_embedding = torch.mean(sequence_embeddings, dim=0)
        return mean_embedding.cpu().numpy()

    def get_batch_embeddings(self, sequences: List[str], batch_size: int = 8, max_length: int = 1024) -> List[List[float]]:
        """Batch embedding with chunking for long sequences."""
        all_embeddings = []
        for i in range(0, len(sequences), batch_size):
            batch_sequences = sequences[i:i + batch_size]
            batch_embeddings = []
            for seq in batch_sequences:
                emb = self.get_embedding(seq, max_length)
                batch_embeddings.append(emb)
            all_embeddings.extend(batch_embeddings)
            if self.device.type in ['cuda', 'mps']:
                torch.cuda.empty_cache() if self.device.type == 'cuda' else None
        return all_embeddings

class ProtBertEmbedder(ProteinEmbedder):
    """Optimized ProtBERT embedder with batch processing."""

    def __init__(self):
        self.device = get_device()
        logger.info(f"Initializing ProtBertEmbedder using device: {self.device}")

        try:
            from transformers import BertModel, BertTokenizer
            self.tokenizer = BertTokenizer.from_pretrained(
                "Rostlab/prot_bert",
                do_lower_case=False
            )
            self.model = BertModel.from_pretrained("Rostlab/prot_bert")
            self.model = self.model.to(self.device)
            self.model.eval()
        except Exception as e:
            logger.error(f"Failed to initialize BERT model: {str(e)}")
            raise EmbeddingError(f"Failed to initialize BERT model: {str(e)}")

    @property
    def vector_size(self) -> int:
        return 1024

    def get_embedding(self, sequence: str) -> List[float]:
        """Single sequence embedding."""
        return self.get_batch_embeddings([sequence], batch_size=1)[0]

    def get_batch_embeddings(self, sequences: List[str], batch_size: int = 8) -> List[List[float]]:
        """Process multiple sequences in batches."""
        all_embeddings = []
        
        for i in range(0, len(sequences), batch_size):
            batch_sequences = sequences[i:i + batch_size]
            
            # Preprocess sequences
            processed_sequences = []
            for seq in batch_sequences:
                seq = " ".join(re.sub(r"[UZOB]", "X", seq))
                processed_sequences.append(seq)
            
            try:
                # Tokenize batch
                inputs = self.tokenizer(
                    processed_sequences,
                    return_tensors='pt',
                    padding=True,
                    truncation=True,
                    max_length=1024
                )
                inputs = {k: v.to(self.device) for k, v in inputs.items()}

                # Generate embeddings
                with torch.no_grad():
                    outputs = self.model(**inputs)
                    embeddings = outputs.last_hidden_state.mean(dim=1)
                    embeddings = embeddings.cpu().numpy()

                # Process each embedding in batch
                for embedding in embeddings:
                    normalized_embedding = normalize_l2(embedding)
                    validated_embedding = self.validate_embedding(normalized_embedding.tolist())
                    all_embeddings.append(validated_embedding)

                # Clear GPU cache
                if self.device.type in ['cuda', 'mps']:
                    torch.cuda.empty_cache() if self.device.type == 'cuda' else None

            except Exception as e:
                logger.error(f"Failed to process batch: {str(e)}")
                raise EmbeddingError(f"Failed to process batch: {str(e)}")
        
        return all_embeddings

class ProtT5Embedder(ProteinEmbedder):
    """Protein embedder using the ProtT5 model."""

    def __init__(self, model_name="Rostlab/prot_t5_xl_half_uniref50-enc"):
        self.device = get_device()
        logger.info(f"Initializing ProtT5Embedder using device: {self.device}")

        try:
            from transformers import T5EncoderModel, T5Tokenizer
            self.tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)
            self.model = T5EncoderModel.from_pretrained(model_name)
            self.model = self.model.to(self.device)
            self.model.eval()
            self._vector_size = self.model.config.d_model
        except Exception as e:
            logger.error(f"Failed to initialize ProtT5 model: {str(e)}")
            raise EmbeddingError(f"Failed to initialize ProtT5 model: {str(e)}")

    @property
    def vector_size(self) -> int:
        return self._vector_size

    def get_embedding(self, sequence: str, max_length: int = 1024) -> List[float]:
        seq = re.sub(r'[^A-Z]', '', sequence.upper())
        seq = re.sub(r'[UZOB]', 'X', seq)
        if len(seq) > max_length - 2:
            seq = seq[:max_length - 2]
        seq = ' '.join(list(seq))
        inputs = self.tokenizer(
            seq,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=max_length,
            add_special_tokens=True
        )
        inputs = {key: val.to(self.device) for key, val in inputs.items()}
        with torch.no_grad():
            outputs = self.model(**inputs)
        last_hidden_state = outputs.last_hidden_state
        attention_mask = inputs['attention_mask']
        valid_positions = attention_mask[0].bool()
        sequence_embeddings = last_hidden_state[0][valid_positions]
        mean_embedding = torch.mean(sequence_embeddings, dim=0)
        embedding_np = mean_embedding.cpu().numpy()
        normalized_embedding = normalize_l2(embedding_np)
        validated_embedding = self.validate_embedding(normalized_embedding.tolist())
        return validated_embedding

    def get_batch_embeddings(self, sequences: List[str], batch_size: int = 8, max_length: int = 1024) -> List[List[float]]:
        all_embeddings = []
        for i in range(0, len(sequences), batch_size):
            batch_sequences = sequences[i:i + batch_size]
            processed_sequences = []
            for seq in batch_sequences:
                seq = re.sub(r'[^A-Z]', '', seq.upper())
                seq = re.sub(r'[UZOB]', 'X', seq)
                if len(seq) > max_length - 2:
                    seq = seq[:max_length - 2]
                processed_sequences.append(' '.join(list(seq)))
            inputs = self.tokenizer(
                processed_sequences,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=max_length,
                add_special_tokens=True
            )
            inputs = {key: val.to(self.device) for key, val in inputs.items()}
            with torch.no_grad():
                outputs = self.model(**inputs)
            last_hidden_state = outputs.last_hidden_state
            attention_mask = inputs['attention_mask']
            for i in range(last_hidden_state.size(0)):
                valid_positions = attention_mask[i].bool()
                sequence_embeddings = last_hidden_state[i][valid_positions]
                mean_embedding = torch.mean(sequence_embeddings, dim=0)
                embedding_np = mean_embedding.cpu().numpy()
                normalized_embedding = normalize_l2(embedding_np)
                validated_embedding = self.validate_embedding(normalized_embedding.tolist())
                all_embeddings.append(validated_embedding)
            if self.device.type in ['cuda', 'mps']:
                torch.cuda.empty_cache() if self.device.type == 'cuda' else None
        return all_embeddings

def get_embedder(model_name: str) -> ProteinEmbedder:
    """Factory function to get the appropriate embedder."""
    try:
        if model_name.lower() == "prot_bert":
            return ProtBertEmbedder()
        elif model_name.lower() == "esm2_small":
            return ESM2Embedder(model_name="facebook/esm2_t12_35M_UR50D")
        elif model_name.lower() == "esm2_large":
            return ESM2Embedder(model_name="facebook/esm2_t36_3B_UR50D")
        elif model_name.lower() == "prot_t5":
            return ProtT5Embedder()
        else:
            raise ValueError(f"Unknown model name: {model_name}")
    except Exception as e:
        logger.error(f"Error creating embedder '{model_name}': {str(e)}")
        raise

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Compute pairwise distances of protein sequences using PLM embeddings (Optimized).')
    parser.add_argument('--input', '-i', required=True, help='Input FASTA file')
    parser.add_argument('--output', '-o', required=False, help='Output CSV file path (optional)')
    parser.add_argument('--model', '-m', default='esm2_small', choices=['prot_bert', 'prot_t5', 'esm2_small', 'esm2_large'], help='Model to use for embeddings')
    parser.add_argument('--batch_size', '-b', type=int, default=8, help='Batch size for processing (default: 8)')
    parser.add_argument('--max_length', type=int, default=1024, help='Maximum sequence length (default: 1024)')
    args = parser.parse_args()

    fasta_file = args.input
    output_file = args.output if args.output else f'{args.model}_dis_matrix_optimized.csv'
    model = args.model
    batch_size = args.batch_size
    max_length = args.max_length

    print(f"Using model: {model}")
    print(f"Batch size: {batch_size}")
    print(f"Max sequence length: {max_length}")
    print(f"Device: {get_device()}")

    embedder = get_embedder(model_name=model)

    # Load all sequences first
    sequences = []
    ids = []
    print("Loading sequences...")
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequences.append(str(record.seq))
        ids.append(str(record.id))
    
    total_sequences = len(sequences)
    print(f"Loaded {total_sequences} sequences")

    # Generate embeddings in batches with progress tracking
    print("Generating embeddings...")
    start_time = time.time()
    
    if hasattr(embedder, 'get_batch_embeddings'):
        # Use optimized batch processing
        all_embeddings = []
        
        for i in range(0, total_sequences, batch_size):
            batch_sequences = sequences[i:i + batch_size]
            batch_ids = ids[i:i + batch_size]
            
            print(f"Processing batch {i//batch_size + 1}/{(total_sequences + batch_size - 1)//batch_size} "
                  f"(sequences {i+1}-{min(i+batch_size, total_sequences)}/{total_sequences})")
            
            try:
                batch_embeddings = embedder.get_batch_embeddings(batch_sequences, batch_size, max_length)
                all_embeddings.extend(batch_embeddings)
                
            except Exception as e:
                print(f"Error processing batch starting at sequence {i+1}: {e}")
                # Process sequences individually as fallback
                for j, seq in enumerate(batch_sequences):
                    try:
                        embedding = embedder.get_embedding(seq)
                        all_embeddings.append(embedding)
                        print(f"  Individual processing: {ids[i+j]}")
                    except Exception as seq_error:
                        print(f"  Failed to process sequence {ids[i+j]}: {seq_error}")
                        # Add zero vector as placeholder
                        all_embeddings.append([0.0] * embedder.vector_size)
    else:
        # Fallback to individual processing
        all_embeddings = []
        for i, (sequence, seq_id) in enumerate(zip(sequences, ids)):
            try:
                print(f'Embedding for Id: {seq_id} ({i+1}/{total_sequences})')
                embedding = embedder.get_embedding(sequence)
                all_embeddings.append(embedding)
            except Exception as e:
                print(f"Error processing {seq_id}: {e}")
                all_embeddings.append([0.0] * embedder.vector_size)

    embedding_time = time.time() - start_time
    print(f"Embedding generation completed in {embedding_time:.2f} seconds")
    print(f"Average time per sequence: {embedding_time/total_sequences:.3f} seconds")

    # Calculate distance matrix with optimizations
    print("Calculating distance matrix...")
    start_time = time.time()
    
    vectors_array = np.array(all_embeddings, dtype=np.float32)  # Use float32 to save memory
    n = len(vectors_array)
    
    # Use symmetric matrix optimization - only calculate upper triangle
    dis_matrix = np.zeros((n, n), dtype=np.float32)
    
    completed_pairs = 0
    total_pairs = n * (n - 1) // 2  # Only unique pairs
    
    for i in range(n):
        for j in range(i + 1, n):  # Only calculate upper triangle
            # Cosine distance calculation
            dot_product = np.dot(vectors_array[i], vectors_array[j])
            norm_i = np.linalg.norm(vectors_array[i])
            norm_j = np.linalg.norm(vectors_array[j])
            
            cosine_similarity = dot_product / (norm_i * norm_j) if (norm_i * norm_j) > 0 else 0
            cosine_distance = 1 - cosine_similarity
            
            dis_matrix[i, j] = cosine_distance
            dis_matrix[j, i] = cosine_distance  # Mirror to lower triangle
            
            completed_pairs += 1
            
            # Show progress every 5% or at completion
            if completed_pairs % max(1, total_pairs // 20) == 0 or completed_pairs == total_pairs:
                progress = (completed_pairs / total_pairs) * 100
                print(f"Distance calculation: {progress:.1f}% ({completed_pairs:,}/{total_pairs:,} pairs)", end='\r')
    
    print()  # New line after progress
    distance_time = time.time() - start_time
    print(f"Distance calculation completed in {distance_time:.2f} seconds")

    # Save results
    print("Saving results...")
    df = pd.DataFrame(dis_matrix, index=ids, columns=ids)
    df.to_csv(output_file, index_label='ID')
    
    total_time = embedding_time + distance_time
    print(f"\nResults saved to: {output_file}")
    print(f"Total runtime: {total_time:.2f} seconds")
    print(f"  - Embedding: {embedding_time:.2f}s ({embedding_time/total_time*100:.1f}%)")
    print(f"  - Distance calculation: {distance_time:.2f}s ({distance_time/total_time*100:.1f}%)")
    print("Done!")