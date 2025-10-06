import logging
import sys
import argparse
import pandas as pd
from Bio import SeqIO
from pymilvus import MilvusClient
from embedders import get_embedder
import concurrent.futures
import threading
import time
from typing import List, Tuple

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('annotation.log')
    ]
)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(
        description='Annotate proteins using vector database'
    )
    parser.add_argument('--input_faa', type=str, required=True, 
                       help='Path to input FAA file to annotate')
    parser.add_argument('--collection', type=str, required=True, 
                       help='Name of the collection to search against')
    parser.add_argument('--output_file', type=str, required=True, 
                       help='Path to output TSV file')
    parser.add_argument('--threshold', type=float, default=0.98, 
                       help='Similarity threshold for annotations')
    parser.add_argument('--model_name', type=str, default="esm2",
                       choices=["prot_bert", "esm2", "esm2_small", "esm2_large", "openai", "prot_t5"],
                       help='Protein embedding model to use')
    parser.add_argument('--embedding_batch_size', type=int, default=16, 
                       help='Batch size for embedding generation')
    parser.add_argument('--num_workers', type=int, default=4, 
                       help='Number of parallel workers for processing')
    return parser.parse_args()

class ProteinAnnotator:
    def __init__(self, args):
        self.args = args
        self.embedder = get_embedder(self.args.model_name)
        self.milvas_client = self._initialize_milvas()
        self.results = []
        self.stats = {
            'total': 0,
            'success': 0,
            'failed': 0,
            'empty': 0,
            'below_threshold': 0
        }
        self._lock = threading.Lock()

    def _initialize_milvas(self):
        try:
            client = MilvusClient('./vector.db')
            # 7. Load the collection
            client.load_collection(
                collection_name=self.args.collection
            )
            return client
        except Exception as e:
            raise ConnectionError(f"Failed to connect to Qdrant database: {str(e)}")

    def process_sequence(self, seq_record):
        """Process a single sequence and find its annotation."""
        seq_id = seq_record.id
        logger.info(f"Processing sequence {seq_id}...")
        seq = str(seq_record.seq).strip('*')
        self.stats['total'] += 1

        if not seq:
            self.stats['empty'] += 1
            return {
                'Query_ID': seq_id,
                'Annotation': 'N/A',
                'Similarity_Score': 0.0,
                'Status': 'empty_sequence'
            }

        try:
            # Get embedding
            embedding = self.embedder.get_embedding(seq)
            
            # Search in database
            res = self.milvas_client.search(
                collection_name=self.args.collection,
                anns_field="vector",
                data=[embedding],
                limit=10,
                search_params={"metric_type": "IP"},
                output_fields=["protein_info"]
            )
            self.stats['success'] += 1
            search_result = [
                    {
                        'Query': seq_record.description,
                        'Annotation': hit['entity']['protein_info'],
                        'Similarity_Score': hit['distance'],
                        'Status': 'success'
                    }
                    for hits in res for hit in hits
            ]
            return search_result
                    

            # if search_result and search_result['distance'] >= self.args.threshold:
            #     self.stats['success'] += 1
            #     return {
            #         'Query_ID': seq_id,
            #         'Annotation': search_result['entity']['protein_info'],
            #         'Similarity_Score': search_result['distance'],
            #         'Status': 'success'
            #     }
            # else:
            #     self.stats['below_threshold'] += 1
            #     return {
            #         'Query_ID': seq_id,
            #         'Annotation':  search_result['entity']['protein_info'] if search_result else 'hypothetical protein',
            #         'Similarity_Score': search_result['distance'] if search_result else 0.0,
            #         'Status': 'below_threshold'
            #     }

        except Exception as e:
            self.stats['failed'] += 1
            logger.error(f"Error processing sequence {seq_id}: {str(e)}")
            return {
                'Query_ID': seq_id,
                'Annotation': 'hypothetical protein',
                'Similarity_Score': 0.0,
                'Status': 'error'
            }

    def _process_sequence_batch(self, sequences_batch: List[Tuple[str, str, str]]) -> List[dict]:
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
                with concurrent.futures.ThreadPoolExecutor(max_workers=self.args.num_workers) as executor:
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
            
            # Process search results for each embedding
            for i, (seq_id, seq, desc) in enumerate(sequences_batch):
                with self._lock:
                    self.stats['total'] += 1
                    
                if not seq:
                    with self._lock:
                        self.stats['empty'] += 1
                    results.append({
                        'Query': desc,
                        'Annotation': 'N/A',
                        'Similarity_Score': 0.0,
                        'Status': 'empty_sequence'
                    })
                    continue
                
                if embeddings[i] is not None:
                    try:
                        # Search in database
                        res = self.milvas_client.search(
                            collection_name=self.args.collection,
                            anns_field="vector",
                            data=[embeddings[i]],
                            limit=10,
                            search_params={"metric_type": "IP"},
                            output_fields=["protein_info"]
                        )
                        
                        batch_results = [
                            {
                                'Query': desc,
                                'Annotation': hit['entity']['protein_info'],
                                'Similarity_Score': hit['distance'],
                                'Status': 'success'
                            }
                            for hits in res for hit in hits
                        ]
                        
                        if batch_results:
                            results.extend(batch_results)
                            with self._lock:
                                self.stats['success'] += 1
                        else:
                            results.append({
                                'Query': desc,
                                'Annotation': 'hypothetical protein',
                                'Similarity_Score': 0.0,
                                'Status': 'no_results'
                            })
                            
                    except Exception as e:
                        logger.error(f"Database search failed for sequence {seq_id}: {str(e)}")
                        with self._lock:
                            self.stats['failed'] += 1
                        results.append({
                            'Query': desc,
                            'Annotation': 'hypothetical protein',
                            'Similarity_Score': 0.0,
                            'Status': 'error'
                        })
                else:
                    with self._lock:
                        self.stats['failed'] += 1
                    results.append({
                        'Query': desc,
                        'Annotation': 'hypothetical protein',
                        'Similarity_Score': 0.0,
                        'Status': 'embedding_failed'
                    })
                        
        except Exception as e:
            logger.error(f"Batch processing failed: {str(e)}")
            with self._lock:
                self.stats['failed'] += len(sequences_batch)
        
        return results

    def run_optimized(self):
        """Execute the optimized annotation pipeline with batch processing."""
        start_time = time.time()
        logger.info("Starting optimized annotation pipeline...")
        
        try:
            # Load all sequences first (single pass)
            sequences = []
            for seq_record in SeqIO.parse(self.args.input_faa, "fasta"):
                seq = str(seq_record.seq).strip('*')
                sequences.append((seq_record.id, seq, seq_record.description))
            
            logger.info(f"Processing {len(sequences)} sequences from {self.args.input_faa}")
            
            # Process in batches
            for batch_idx, i in enumerate(range(0, len(sequences), self.args.embedding_batch_size)):
                batch = sequences[i:i + self.args.embedding_batch_size]
                
                # Log start of batch processing
                start_idx = i + 1
                end_idx = min(i + self.args.embedding_batch_size, len(sequences))
                logger.info(f'Starting annotation batch {batch_idx + 1}: sequences {start_idx}-{end_idx}')
                
                batch_start_time = time.time()
                batch_results = self._process_sequence_batch(batch)
                batch_time = time.time() - batch_start_time
                
                logger.info(f'Completed annotation batch {batch_idx + 1} in {batch_time:.1f}s ({len(batch_results)} results)')
                
                # Add to results
                self.results.extend(batch_results)
                
                # Progress logging with time estimation
                progress_percentage = (self.stats['total'] / len(sequences)) * 100
                elapsed_time = time.time() - start_time
                
                if self.stats['total'] > 0:
                    avg_time_per_seq = elapsed_time / self.stats['total']
                    remaining_sequences = len(sequences) - self.stats['total']
                    estimated_remaining_time = avg_time_per_seq * remaining_sequences
                    
                    # Format time estimates
                    elapsed_str = f"{elapsed_time:.1f}s"
                    remaining_str = f"{estimated_remaining_time:.1f}s"
                    
                    logger.info(f'Annotation Progress: {self.stats["total"]}/{len(sequences)} sequences ({progress_percentage:.1f}%) | Elapsed: {elapsed_str} | ETA: {remaining_str}')
                else:
                    logger.info(f'Annotation Progress: {self.stats["total"]}/{len(sequences)} sequences ({progress_percentage:.1f}%)')

            # Save results
            df = pd.DataFrame(self.results)
            df.to_csv(self.args.output_file, sep='\t', index=False)
            logger.info(f"Results saved to {self.args.output_file}")

            # Final summary
            total_time = time.time() - start_time
            avg_time_per_seq = total_time / len(sequences) if len(sequences) > 0 else 0
            logger.info(f"Annotation Complete: {self.stats['total']}/{len(sequences)} sequences processed (100.0%)")
            logger.info(f"Total processing time: {total_time:.1f}s | Average: {avg_time_per_seq:.3f}s per sequence")

            # Print summary
            self._print_summary()

        except Exception as e:
            logger.error(f"Pipeline failed: {str(e)}")
            raise

    def run(self):
        """Execute the annotation pipeline."""
        logger.info("Starting annotation pipeline...")
        
        try:
            # Process sequences one by one
            for seq_record in SeqIO.parse(self.args.input_faa, "fasta"):
                rets = self.process_sequence(seq_record)
                # self.results.extend(rets)
                if isinstance(rets, list):
                    self.results.extend(rets)
                else:
                    self.results.append(rets)
                # print(self.results)            
                # Log progress every 100 sequences
                if len(self.results) % 100 == 0:
                    logger.info(f"Processed a batch ...")

            # Save results
            df = pd.DataFrame(self.results)
            df.to_csv(self.args.output_file, sep='\t', index=False)
            logger.info(f"Results saved to {self.args.output_file}")

            # Print summary
            self._print_summary()

        except Exception as e:
            logger.error(f"Pipeline failed: {str(e)}")
            raise

    def _print_summary(self):
        """Print summary statistics of the annotation results."""
        logger.info("\nAnnotation Summary:")
        logger.info(f"Total sequences processed: {self.stats['total']}")
        logger.info(f"Successful annotations: {self.stats['success']}")
        logger.info(f"Below threshold: {self.stats['below_threshold']}")
        logger.info(f"Failed sequences: {self.stats['failed']}")
        logger.info(f"Empty sequences: {self.stats['empty']}")

def main():
    args = parse_args()
    logger.info(f"Loading embedding model: {args.model_name}")
    annotator = ProteinAnnotator(args)
    logger.info("Annotator initialized successfully")
    annotator.run_optimized()

if __name__ == '__main__':
    main()
