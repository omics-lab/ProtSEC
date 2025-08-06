import argparse
from Bio import SeqIO
import pickle
from embedder import ProteinEmbedder
import numpy as np
from scipy.fft import ifft
import time


def parse_args():
    parser = argparse.ArgumentParser(
        description='Annotate proteins using vector database'
    )
    parser.add_argument('--input_faa', type=str, required=True, help='Path to input FAA file to annotate')
    parser.add_argument('--dim_reduct', type=str, default='MDS', choices=['UMAP', 't-SNE', 'MDS'],
                        help='Algorithm for dimensionality reduction (default: t-SNE)')
    parser.add_argument('--dist_func', type=str, default='ASMP', choices=['SMS', 'ASMP', 'SNN'], 
                        help='Distance function for computing distance matrix (default: SMS)')
    parser.add_argument('--top_hit', type=int, default=1, help='Number of top hits to return')
    parser.add_argument('--db', type=str, required=True, help='path to the precomputed database')
    parser.add_argument('--out', type=str, required=True, help='Path to output TSV file')

    return parser.parse_args()

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


def calculate_score(F1, F2):
    # Compute the normalized cross-power spectrum
    cross_power_spectrum = F1 * np.conj(F2)
    cross_power_spectrum /= (np.abs(cross_power_spectrum))

    # Compute inverse FFT to obtain the correlation
    correlation = ifft(cross_power_spectrum)

    # Find the peak value and index
    peak_value = np.max(np.abs(correlation))
    # shift_index = np.argmax(np.abs(correlation))
    
    return peak_value





def main():
    args = parse_args()
    
    # Open the file in binary read mode ('rb')
    with open(args.db, 'rb') as f:
        # Load the data from the file
        list_of_vectors = pickle.load(f)
    
    # print(list_of_vectors)
    encoder = ProteinEmbedder(args.dim_reduct, args.dist_func)
    print(f"Using dimensionality reduction method: {args.dim_reduct}")
    print(f"Using distance function: {args.dist_func}")
    
    print(f"Using {args.top_hit} top hits")
    sequences = list(SeqIO.parse(args.input_faa, 'fasta'))
    print(f"Found {len(sequences)} sequences to annotate")

    # Clear output file
    with open(args.out, 'w') as f:
        pass
    
    
    for seq_record in sequences:
        start_time = time.time()
        res = encode_sequence(seq_record, encoder)
        if res is None:
            continue
        # Calculate the score of the sequence with all the vectors
        scores = [calculate_score(res['v'], v['v']) for v in list_of_vectors]
        # scores = sorted(scores, reverse=True)


        scores_copy = np.array(scores, copy=True)
        
        top_indices = []
        # Find top k elements
        for _ in range(args.top_hit):
            idx = np.argmax(scores_copy)
            top_indices.append(idx)
            scores_copy[idx] = float('-inf')
        
        # Write the results to a TSV file
        with open(args.out, 'a') as f:
            for idx in top_indices:
                f.write(f"{seq_record.description}\t{list_of_vectors[idx]['info']}\t{scores[idx]}\n")
        
        elapsed_time = time.time() - start_time

        print(f"Annotated {seq_record.id} Time elapsed: {elapsed_time:.4f} seconds")

        
    print(f"Results written to {args.out}")

    

if __name__ == '__main__':
    main()

# python3 annotate.py \
#     --input_faa ./data/QUERY.fasta \
#     --db ./DB/mmseq2_db.pkl \
#     --out ./data/mmseq2_result.tsv