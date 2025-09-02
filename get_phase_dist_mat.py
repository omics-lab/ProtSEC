from embedder import ProteinEmbedder
from Bio import SeqIO

import numpy as np
from scipy.fft import ifft
import argparse
import pandas as pd


def calculate_distance(F1, F2):
    # Compute the normalized cross-power spectrum
    cross_power_spectrum = F1 * np.conj(F2)
    cross_power_spectrum /= (np.abs(cross_power_spectrum))

    # Compute inverse FFT to obtain the correlation
    correlation = ifft(cross_power_spectrum)

    # Find the peak value (correlation score)
    peak_value = np.max(np.abs(correlation))
    
    # Convert correlation to distance (1 - correlation)
    distance = 1 - peak_value
    
    return distance


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Compute pairwise distances of protein sequences from phase correlation')
    parser.add_argument('--input', '-i', required=True, help='Input FASTA file')
    parser.add_argument('--output', '-o', required=True, help='Output CSV file path')
    parser.add_argument('--dim', '-n', type=int, default=1024, help='FFT dimension (default: 1024)')
    args = parser.parse_args()


    dataset_path = args.input
    output_csv = args.output
    fft_dim = args.dim
    embedder = ProteinEmbedder(n=fft_dim)

    ids = []
    features = []

    print("Encoding sequences...")
    for i, record in enumerate(SeqIO.parse(dataset_path, 'fasta'), 1):
        v = embedder.encode(str(record.seq))
        ids.append(str(record.id))
        features.append(v)
        
        # Show progress every 10 sequences or for the last sequence
        if i % 10 == 0 or i == 1:
            print(f"Encoded {i} sequences...", end='\r')
    
    print(f"\nCompleted encoding {len(ids)} sequences")

    n = len(ids)
    print(f'Number of sequences: {n}')

    distances = np.zeros((n, n))

    print("Calculating pairwise distances...")
    total_pairs = n * n
    completed_pairs = 0

    for i in range(n):
        for j in range(n):
            distance = calculate_distance(features[i], features[j])
            distances[i, j] = round(distance, 5)
            completed_pairs += 1
            
            # Show progress every 100 pairs or for significant milestones
            if completed_pairs % 100 == 0 or completed_pairs == total_pairs:
                progress = (completed_pairs / total_pairs) * 100
                print(f"Progress: {completed_pairs}/{total_pairs} pairs ({progress:.1f}%)", end='\r')
    
    print(f"\nCompleted calculating all {total_pairs} pairwise distances")

    # Step 4: Create DataFrame and save as CSV
    print("Saving results to CSV...")
    df = pd.DataFrame(distances, index=ids, columns=ids)
    df.to_csv(output_csv, index_label='ID')
    print(f"Results saved to: {output_csv}")

# python3 get_phase_dist_mat.py -n 512 -i phosphatase.fa -o score_matrix.csv
