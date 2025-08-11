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
        # Set device explicitly
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

    def validate_embedding(self, embedding: List[float]) -> List[float]:
        """Validate embeddings to ensure they are well-formed."""
        if not embedding:
            raise EmbeddingError("Generated embedding is empty")
        if not all(isinstance(x, float) for x in embedding):
            raise EmbeddingError("Non-float values in embedding")
        if any(np.isnan(x) or np.isinf(x) for x in embedding):
            raise EmbeddingError("NaN or Inf values in embedding")
        return embedding


class ProtBertEmbedder(ProteinEmbedder):
    """Protein embedder using the ProtBERT model."""

    def __init__(self):
        # Ensure device is set before model initialization
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
        try:
            # Preprocess sequence: replace rare amino acids with X and space out
            sequence = " ".join(re.sub(r"[UZOB]", "X", sequence))

            # Tokenize and encode
            encoded_input = self.tokenizer(sequence, return_tensors='pt')
            encoded_input = {k: v.to(self.device) for k, v in encoded_input.items()}

            # Generate embeddings
            with torch.no_grad():
                outputs = self.model(**encoded_input)
                embeddings = outputs.last_hidden_state.mean(dim=1)
                embeddings = embeddings.cpu()

            # Normalize and convert to list
            normalized_embedding = normalize_l2(embeddings.numpy()[0])
            embedding_list = normalized_embedding.tolist()

            return self.validate_embedding(embedding_list)

        except Exception as e:
            logger.error(f"Failed to generate embedding: {str(e)}")
            raise EmbeddingError(f"Failed to generate embedding: {str(e)}")


class ProtT5Embedder(ProteinEmbedder):
    """Protein embedder using the ProtT5 model."""
    def __init__(self, model_name="Rostlab/prot_t5_xl_half_uniref50-enc"):
        # Ensure device is set before model initialization
        self.device = get_device()
        logger.info(f"Initializing ProtT5Embedder using device: {self.device}")

        try:
            from transformers import T5EncoderModel, T5Tokenizer
            self.tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)
            self.model = T5EncoderModel.from_pretrained(model_name)
            self.model = self.model.to(self.device)
            self.model = self.model.eval()
        except Exception as e:
            logger.error(f"Failed to initialize T5 model: {str(e)}")
            raise EmbeddingError(f"Failed to initialize T5 model: {str(e)}")
    @property
    def vector_size(self) -> int:
        return 1024

    def get_embedding(self, sequence: str) -> List[float]:
        try:
            # Preprocess sequence: replace rare amino acids with X and space out
            sequence = " ".join(re.sub(r"[UZOB]", "X", sequence))

            # Tokenize sequence
            token_encoding = self.tokenizer.encode_plus(sequence, add_special_tokens=True, return_tensors='pt')
            input_ids = token_encoding['input_ids'].to(self.device)
            attention_mask = token_encoding['attention_mask'].to(self.device)

            # Generate embeddings
            with torch.no_grad():
                embedding_repr = self.model(input_ids, attention_mask=attention_mask)

            # Extract embeddings (removing special tokens)
            emb = embedding_repr.last_hidden_state[0, :len(sequence)]
            result = emb.mean(dim=0).detach().cpu().numpy()

            # Normalize and convert to list
            normalized_embedding = normalize_l2(result)
            embedding_list = normalized_embedding.tolist()

            return self.validate_embedding(embedding_list)

        except Exception as e:
            logger.error(f"Failed to generate embedding: {str(e)}")
            raise EmbeddingError(f"Failed to generate embedding: {str(e)}")


class ESM2Embedder(ProteinEmbedder):
    """Protein embedder using the ESM2 model."""

    def __init__(self, model_name="facebook/esm2_t12_35M_UR50D"):
        # IMPORTANT: Set device explicitly before any other initialization
        self.device = get_device()
        logger.info(f"Initializing ESM2Embedder using device: {self.device}")

        try:
            from transformers import AutoTokenizer, AutoModel
            self.model_name = model_name
            self.tokenizer = AutoTokenizer.from_pretrained(model_name)
            self.model = AutoModel.from_pretrained(model_name)
            self.model = self.model.to(self.device)
            self.model.eval()

            # Get embedding size from model config
            self._vector_size = self.model.config.hidden_size

        except Exception as e:
            logger.error(f"Failed to initialize ESM2 model: {str(e)}")
            raise EmbeddingError(f"Failed to initialize ESM2 model: {str(e)}")

    @property
    def vector_size(self) -> int:
        return self._vector_size

    def _get_chunk_embedding(self, sequence_chunk: str, max_length: int = 1024) -> np.ndarray:
        """Process a single chunk of sequence and return its embedding."""
        # Tokenize and encode with attention mask to handle padding properly
        inputs = self.tokenizer(
            sequence_chunk,
            return_tensors="pt",
            padding="max_length",
            truncation=True,
            max_length=max_length,
            add_special_tokens=True
        )

        # Move to device
        # print(f"DEBUG: self.device = {self.device}")
        inputs = {key: val.to(self.device) for key, val in inputs.items()}

        # Generate embeddings
        with torch.no_grad():
            outputs = self.model(**inputs)

        # Get sequence embeddings
        # For ESM2, we can use the last hidden state
        last_hidden_state = outputs.last_hidden_state

        # Mean pooling - take average of all token embeddings
        emb = torch.mean(last_hidden_state, dim=1)

        return emb.cpu().numpy()[0]

    def get_embedding(self, sequence: str, max_length: int = 1024) -> List[float]:
        try:
            # Preprocess sequence - some ESM models require specific formatting
            # Remove any whitespace and non-amino acid characters
            sequence = re.sub(r'[^A-Z]', '', sequence.upper())

            # Replace rare amino acids with X
            sequence = re.sub(r'[UZOB]', 'X', sequence)

            # Calculate how many chunks we need (subtract 2 to account for special tokens [CLS] and [SEP])
            effective_max_length = max_length - 2

            # If sequence is shorter than max_length, process it directly
            if len(sequence) <= effective_max_length:
                # print('hi...')
                embedding = self._get_chunk_embedding(sequence, max_length)
                normalized_embedding = normalize_l2(embedding)
                return self.validate_embedding(normalized_embedding.tolist())

            # For long sequences, split into chunks and process each separately
            embeddings = []

            # Process complete chunks of size effective_max_length
            for i in range(0, len(sequence), effective_max_length):
                chunk = sequence[i:i + effective_max_length]
                if len(chunk) > 0:  # Only process non-empty chunks
                    chunk_embedding = self._get_chunk_embedding(chunk, max_length)
                    embeddings.append(chunk_embedding)

            # Average all chunk embeddings
            if embeddings:
                avg_embedding = np.mean(embeddings, axis=0)
                normalized_embedding = normalize_l2(avg_embedding)
                return self.validate_embedding(normalized_embedding.tolist())
            else:
                raise EmbeddingError("No valid chunks were processed")

        except Exception as e:
            logger.error(f"Failed to generate embedding: {str(e)}")
            raise EmbeddingError(f"Failed to generate embedding: {str(e)}")


class OpenAIEmbedder(ProteinEmbedder):
    """Protein embedder using the OpenAI GPT model."""

    def __init__(self, model_name="text-embedding-3-large"):
        time.sleep(0.001)
        try:
            from openai import OpenAI
            self._vector_size = 3072
            self.client = OpenAI()
            self.model_name = model_name

        except Exception as e:
            logger.error(f"Failed to initialize OpenAI clinet: {str(e)}")
            raise EmbeddingError(f"Failed to initialize OpenAI client: {str(e)}")

    @property
    def vector_size(self) -> int:
        return self._vector_size

    def get_embedding(self, sequence: str) -> List[float]:
        try:
            sequence = " ".join(re.sub(r"[UZOB]", "X", sequence))

            response = self.client.embeddings.create(
                input=sequence,
                model=self.model_name
            )
            return response.data[0].embedding

        except Exception as e:
            logger.error(f"Failed to generate embedding: {str(e)}")
            raise EmbeddingError(f"Failed to generate embedding: {str(e)}")


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
        elif model_name.lower() == "openai":
            return OpenAIEmbedder()
        else:
            raise ValueError(f"Unknown model name: {model_name}")
    except Exception as e:
        logger.error(f"Error creating embedder '{model_name}': {str(e)}")
        raise


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Compute pairwise distances of protein sequences using PLM embeddings.')
    parser.add_argument('--input', '-i', required=True, help='Input FASTA file')
    # parser.add_argument('--output', '-o', required=True, help='Output CSV file path')
    parser.add_argument('--model', '-m', default='esm2_small', choices=['prot_bert', 'esm2_small', 'esm2_large', 'prot_t5'], help='Model to use for embeddings')
    args = parser.parse_args()

    fasta_file = args.input
    # output_file = args.output
    model = args.model

    embedder = get_embedder(model_name=model)

    vectors, ids = [], []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequence = str(record.seq)
        print(f'Embedding for Id: {record.id}')
        embedding = embedder.get_embedding(sequence)
        id = str(record.id)
        
        vectors.append(embedding)
        ids.append(id)
    

    dis_matrix = cdist(np.array(vectors), np.array(vectors), metric='cosine')

    # Step 4: Create DataFrame and save as CSV
    df = pd.DataFrame(dis_matrix, index=ids, columns=ids)
    df.to_csv(f'{model}_dis_matrix.csv', index_label='ID')
    print(f'Distance matrix saved to {model}_dis_matrix.csv')
    print("Done!")