import argparse
from Bio import SeqIO
import pickle
from embedder import ProteinEmbedder
import numpy as np
from fastdtw import fastdtw
from scipy.stats import pearsonr
from dtaidistance import dtw_ndim


def parse_args():
    parser = argparse.ArgumentParser(
        description='Annotate proteins using vector database'
    )
    parser.add_argument('--input_faa', type=str, required=True, help='Path to input FAA file to annotate')
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


# def calculate_score(v1, v2):
#     """Calculate the score between two vectors."""
#     inner_product = np.vdot(v1, v2)
#     return np.abs(inner_product)

# def calculate_score(F1, F2):
#     l = max(len(F1), len(F2))
#     F1 = np.pad(F1, (0, l - len(F1)), 'constant', constant_values=0+0j)
#     F2 = np.pad(F2, (0, l - len(F2)), 'constant', constant_values=0+0j)
#     # print(F1)

#     # Compute the normalized cross-power spectrum
#     cross_power_spectrum = F1 * np.conj(F2)
#     cross_power_spectrum /= np.abs(cross_power_spectrum)

#     # Compute inverse FFT to obtain the correlation
#     correlation = np.fft.ifft(cross_power_spectrum)
#     correlation = np.abs(correlation)
#     # Find the peak value and index
#     peak_value = np.max(np.abs(correlation))
#     # shift_index = np.argmax(np.abs(correlation))
    
#     return peak_value

def calculate_score(X, Y):
    """Compute similarity score using normalized cross-spectral function."""
    # Compute FFT of both signals
    # X = np.fft.fft(x)
    # Y = np.fft.fft(y)

    # Compute the cross-spectrum Sxy = X * conj(Y)
    Sxy = X * np.conj(Y)

    # Compute auto power spectral densities (PSDs)
    Sxx = np.abs(X)
    Syy = np.abs(Y)

    # Compute magnitude coherence
    Cxy = np.abs(Sxy) / np.sqrt(Sxx * Syy)

    # Return the mean similarity score across all frequencies
    return np.mean(Cxy)

# def calculate_score(complex_vector1, complex_vector2):
#     v1 = np.column_stack((np.real(complex_vector1), np.imag(complex_vector1)))
#     v2 = np.column_stack((np.real(complex_vector2), np.imag(complex_vector2)))
#     d = dtw_ndim.distance(v1, v2)
#     return d

# def calculate_score(complex_vector1, complex_vector2):
#     """Calculate the score between two vectors."""
#     max_len = max(len(complex_vector1), len(complex_vector2))
#     padded1 = np.zeros(max_len, dtype=complex)
#     padded2 = np.zeros(max_len, dtype=complex)
    
#     padded1[:len(complex_vector1)] = complex_vector1
#     padded2[:len(complex_vector2)] = complex_vector2
    
#     # Compute cross-power spectrum
#     cross_power = np.conj(padded1) * padded2
#     cross_power_norm = cross_power / (np.abs(cross_power) + 1e-10)
    
#     # Inverse DFT to get correlation
#     correlation = np.fft.ifft(cross_power_norm)
#     correlation_mag = np.abs(correlation)
    
#     # Find the peak
#     shift = np.argmax(correlation_mag)
#     if shift > max_len // 2:
#         shift = shift - max_len  # Handle negative shift
    
#     similarity = np.max(correlation_mag)
#     return similarity




def main():
    args = parse_args()
    
    # Open the file in binary read mode ('rb')
    with open(args.db, 'rb') as f:
        # Load the data from the file
        list_of_vectors = pickle.load(f)
    
    # print(list_of_vectors)
    encoder = ProteinEmbedder()
    sequences = list(SeqIO.parse(args.input_faa, 'fasta'))
    print(f"Found {len(sequences)} sequences to annotate")
    
    for seq_record in sequences:
        res = encode_sequence(seq_record, encoder)
        if res is None:
            continue
        # Calculate the score of the sequence with all the vectors
        scores = [calculate_score(res['v'], v['v']) for v in list_of_vectors]
        scores_copy = np.array(scores, copy=True)
        # scores = [(calculate_score(res['v'], v['v']), v['info']) for v in list_of_vectors]
        top_indices = []
        # Find top k elements
        for _ in range(5):
            idx = np.argmax(scores_copy)
            top_indices.append(idx)
            scores_copy[idx] = float('-inf')
        
        # for idx in top_indices:
        #     print(f"{seq_record.description}\t{list_of_vectors[idx]['info']}\t{scores[idx][0]}\t{res['length']}") 
        
    # Write the results to a TSV file
        with open(args.out, 'a') as f:
            for idx in top_indices:
                f.write(f"{seq_record.id}\t{list_of_vectors[idx]['info']}\t{scores[idx]}\n")

        
    print(f"Results written to {args.out}")

    

if __name__ == '__main__':
    main()

#  python3 annotate.py --input_faa ./example_data/ATP-dependent_DNA_helicase_RecG_only.fasta  --db ./DB/protein_data.pkl --out atp.tsv
