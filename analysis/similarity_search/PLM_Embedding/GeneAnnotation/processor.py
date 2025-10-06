import logging
from dataclasses import dataclass
from typing import List, Tuple
from Bio import SeqIO
import concurrent.futures
import threading
from queue import Queue
import time

logger = logging.getLogger(__name__)


@dataclass
class ProcessingStats:
    """Statistics for tracking progress and errors."""
    total_sequences: int = 0
    processed_sequences: int = 0
    failed_sequences: int = 0
    empty_sequences: int = 0

    def increment_processed(self):
        self.processed_sequences += 1

    def increment_failed(self):
        self.failed_sequences += 1

    def increment_empty(self):
        self.empty_sequences += 1


class ProteinProcessor:
    """Handles processing of protein sequences and insertion into database."""

    def __init__(self, embedder, db, batch_size=50, num_workers=4, embedding_batch_size=16):
        """Initialize the processor with necessary components."""
        self.embedder = embedder
        self.db = db
        self.batch_size = batch_size
        self.embedding_batch_size = embedding_batch_size
        self.num_workers = num_workers
        self.points = []
        self.stats = ProcessingStats()
        self._lock = threading.Lock()
    

    def _process_sequence_batch(self, sequences_batch: List[Tuple[int, str, str]]) -> List[dict]:
        """Process a batch of sequences in parallel."""
        results = []
        
        # Extract sequences for batch embedding
        seqs_only = [seq for _, seq, _ in sequences_batch]
        seq_ids = [seq_id for seq_id, _, _ in sequences_batch]
        descriptions = [desc for _, _, desc in sequences_batch]
        
        try:
            # Get batch embeddings if embedder supports it
            if hasattr(self.embedder, 'get_batch_embeddings'):
                embeddings = self.embedder.get_batch_embeddings(seqs_only)
            else:
                # Fallback to individual embeddings with threading
                embeddings = []
                with concurrent.futures.ThreadPoolExecutor(max_workers=self.num_workers) as executor:
                    future_to_seq = {executor.submit(self.embedder.get_embedding, seq): i 
                                   for i, seq in enumerate(seqs_only)}
                    
                    embeddings = [None] * len(seqs_only)
                    for future in concurrent.futures.as_completed(future_to_seq):
                        seq_idx = future_to_seq[future]
                        try:
                            embeddings[seq_idx] = future.result()
                        except Exception as e:
                            logger.error(f"Failed to process sequence {seq_ids[seq_idx]}: {str(e)}")
                            embeddings[seq_idx] = None
            
            # Create results
            for i, (seq_id, seq, desc) in enumerate(sequences_batch):
                if embeddings[i] is not None:
                    results.append({
                        "id": seq_id,
                        "vector": embeddings[i],
                        "protein_info": desc,
                        "sequence_length": len(seq)
                    })
                    with self._lock:
                        self.stats.increment_processed()
                else:
                    with self._lock:
                        self.stats.increment_failed()
                        
        except Exception as e:
            logger.error(f"Batch processing failed: {str(e)}")
            with self._lock:
                self.stats.failed_sequences += len(sequences_batch)
        
        return results

    def process_fasta_file_optimized(self, fasta_path: str) -> ProcessingStats:
        """Optimized version with batch processing and parallelization."""
        start_time = time.time()
        
        # Load all sequences first (single pass)
        sequences = []
        for i, seq_record in enumerate(SeqIO.parse(fasta_path, "fasta"), 1):
            seq = str(seq_record.seq)
            if not seq:
                self.stats.increment_empty()
                continue
            sequences.append((i, seq, seq_record.description))
        
        self.stats.total_sequences = len(sequences) + self.stats.empty_sequences
        logger.info(f"Processing {len(sequences)} valid sequences from {fasta_path}")
        
        # Process in batches
        for batch_idx, i in enumerate(range(0, len(sequences), self.embedding_batch_size)):
            batch = sequences[i:i + self.embedding_batch_size]
            
            # Log start of batch processing
            start_idx = i + 1
            end_idx = min(i + self.embedding_batch_size, len(sequences))
            logger.info(f'Starting embedding batch {batch_idx + 1}: sequences {start_idx}-{end_idx}')
            
            batch_start_time = time.time()
            batch_results = self._process_sequence_batch(batch)
            batch_time = time.time() - batch_start_time
            
            logger.info(f'Completed embedding batch {batch_idx + 1} in {batch_time:.1f}s ({len(batch_results)} successful)')
            
            # Add to points
            self.points.extend(batch_results)
            
            # Progress logging after each embedding batch with time estimation
            progress_percentage = (self.stats.processed_sequences / len(sequences)) * 100
            elapsed_time = time.time() - start_time
            
            if self.stats.processed_sequences > 0:
                avg_time_per_seq = elapsed_time / self.stats.processed_sequences
                remaining_sequences = len(sequences) - self.stats.processed_sequences
                estimated_remaining_time = avg_time_per_seq * remaining_sequences
                
                # Format time estimates
                elapsed_str = f"{elapsed_time:.1f}s"
                remaining_str = f"{estimated_remaining_time:.1f}s"
                
                logger.info(f'Embedding Progress: {self.stats.processed_sequences}/{len(sequences)} sequences ({progress_percentage:.1f}%) | Elapsed: {elapsed_str} | ETA: {remaining_str}')
            else:
                logger.info(f'Embedding Progress: {self.stats.processed_sequences}/{len(sequences)} sequences ({progress_percentage:.1f}%)')
            
            # Upload if batch is full
            if len(self.points) >= self.batch_size:
                try:
                    self.db.upload_points(self.points)
                    self.points = []
                    logger.info(f'Database Upload: {len(batch_results)} sequences uploaded to database')
                except Exception as e:
                    logger.error(f"Failed to upload points: {str(e)}")
        
        # Upload remaining points
        if self.points:
            try:
                self.db.upload_points(self.points)
                logger.info(f"Final Upload: {len(self.points)} sequences uploaded to database")
            except Exception as e:
                logger.error(f"Failed to upload final batch of points: {str(e)}")
        
        # Final progress summary
        total_time = time.time() - start_time
        avg_time_per_seq = total_time / len(sequences) if len(sequences) > 0 else 0
        logger.info(f"Embedding Complete: {self.stats.processed_sequences}/{len(sequences)} sequences successfully embedded (100.0%)")
        logger.info(f"Total processing time: {total_time:.1f}s | Average: {avg_time_per_seq:.3f}s per sequence")
        
        return self.stats

    def process_fasta_file(self, fasta_path: str) -> ProcessingStats:
        """Process all sequences in a FASTA file and insert into database."""
        # Count total sequences
        self.stats.total_sequences = sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))
        logger.info(f"Processing {self.stats.total_sequences} sequences from {fasta_path}")
        
        for i, seq_record in enumerate(SeqIO.parse(fasta_path, "fasta"), 1):
            seq = str(seq_record.seq)
            if not seq:
                self.stats.increment_empty()
                continue
            try:
                v = self.embedder.get_embedding(seq)
                self.points.append(
                    {
                        "id": i,
                        "vector": v,
                        "protein_info": seq_record.description,
                        "sequence_length": len(seq)
                    }
                )
                self.stats.increment_processed()

            except Exception as e:
                self.stats.increment_failed()
                logger.error(f"Failed to process sequence: {seq_record.id} ... - {str(e)}")
            
            if len(self.points) >= self.batch_size:
                try:
                    self.db.upload_points(self.points)
                    
                    self.points = []
                    logger.info(f'Successfully Processed {self.stats.processed_sequences}/{self.stats.total_sequences}')
                except Exception as e:
                    logger.error(f"Failed to upload points: {str(e)}")
        
        # Upload remaining sequences in the final batch if any
        if self.points:
            try:
                self.db.upload_points(self.points)
                logger.info(f"Uploaded final batch of {len(self.points)} points")
            except Exception as e:
                logger.error(f"Failed to upload final batch of points: {str(e)}")

        return self.stats
