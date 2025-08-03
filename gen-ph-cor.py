from embedder import ProteinEmbedder
from Bio import SeqIO

import numpy as np
from scipy.fft import ifft
import argparse
import pandas as pd


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


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Compute pairwise phase correlation scores for protein sequences.')
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

    for record in SeqIO.parse(dataset_path, 'fasta'):
        v = embedder.encode(str(record.seq))
        ids.append(str(record.id))
        features.append(v)

    n = len(ids)
    print(f'Number of sequences: {n}')

    scores = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            score = calculate_score(features[i], features[j])
            scores[i, j] = round(score, 5)

    # Step 4: Create DataFrame and save as CSV
    df = pd.DataFrame(scores, index=ids, columns=ids)
    df.to_csv(output_csv, index_label='ID')

# python3 gen-ph-cor.py -n 512 -i phosphatase.fa -o score_matrix.csv
