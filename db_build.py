from pathlib import Path
import argparse
import multiprocessing
from Bio import SeqIO
from embedder import ProteinEmbedder
from concurrent.futures import ThreadPoolExecutor, as_completed
import pickle


def parse_args() -> argparse.Namespace:
    """Parse command line arguments with input validation."""
    parser = argparse.ArgumentParser(
        description='Build protein vector database from FASTA file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--fasta_path', type=str, required=True, help='Path to input FASTA file')
    # parser.add_argument('--dim', type=int, required=True, help='Dimensionality of the embeddings')
    parser.add_argument('--num_workers', type=int, default=multiprocessing.cpu_count(), 
                       help='Number of worker threads to use (default: number of CPU cores)')
    parser.add_argument('--dim_reduct', type=str, default='MDS', choices=['UMAP', 't-SNE', 'MDS'],
                        help='Algorithm for dimensionality reduction (default: t-SNE)')
    parser.add_argument('--dist_func', type=str, default='SMS', choices=['SMS', 'ASMP', 'SNN'], 
                        help='Distance function for computing distance matrix (default: SMS)')
    parser.add_argument('--db', type=str, default='protein_data.pkl', help='Filename of the database')
    args = parser.parse_args()

    # Validate arguments
    if not Path(args.fasta_path).exists():
        raise FileNotFoundError(f"FASTA file not found: {args.fasta_path}")
    return args


def encode_sequence(seq_record, protein_encoder):
    """Process a single sequence with the protein encoder."""
    sequence = str(seq_record.seq).strip('*')
    if not sequence:
        return None
        
    vector = protein_encoder.encode(sequence)
    
    return {
        'v': vector,
        'info': seq_record.description,
        'length': len(sequence)
    }


def main():
    args = parse_args()
    fasta_path = args.fasta_path
    # dim = args.dim
    num_workers = args.num_workers

    # dimension reduction method
    dim_reduct = args.dim_reduct
    print(f"Using dimensionality reduction method: {dim_reduct}")

    dist_func = args.dist_func
    print(f"Using distance function: {dist_func}")
    
    print(f"Using {num_workers} worker threads")
    
    # Initialize the encoder
    protein_encoder = ProteinEmbedder(dim_reduct, dist_func)
    
    # Read all sequences from the FASTA file first
    sequences = list(SeqIO.parse(args.fasta_path, "fasta"))
    total_sequences = len(sequences)
    print(f"Found {total_sequences} sequences to encode")
    
    list_of_vectors = []
    cnt = 0
    # Process sequences in parallel using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submit all encoding tasks to the executor
        futures = [
            executor.submit(encode_sequence, seq, protein_encoder) for seq in sequences
        ]
        
        # Process results as they complete
        for future in as_completed(futures):
            try:
                result = future.result()
                if result:
                    list_of_vectors.append(result)
                    cnt += 1
                    if cnt % 100 == 0:  # Print progress every 100 sequences
                        print(f"Processed {cnt}/{total_sequences} sequences")
            except Exception as exc:
                print(f"Seq generated an exception: {exc}")
    
    print(f"Successfully encoded {len(list_of_vectors)} sequences")
    
    db_path = Path('DB')  # Replace 'DB' with your directory path
    if not db_path.exists():
        db_path.mkdir(parents=True, exist_ok=True)
    
    with open(f'DB/{args.db}', 'wb') as f:
        pickle.dump(list_of_vectors, f)
    # print(list_of_vectors)
    print(f"Database saved to: DB/{args.db}")

if __name__ == '__main__':
    main()
