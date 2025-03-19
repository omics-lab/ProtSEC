from scipy.fft import fft
import numpy as np
from utility import aa_to_complex

def normalize(x: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(x)
    return x / (norm if norm > 0 else 1)


class ProteinEmbedder():

    def __init__(self):
        pass

    def encode(self, sequence):
        """
        Encode a protein sequence into a complex vector.

        Args:
            sequence: Protein sequence as string

        Returns:
           complex vector representing the sequence
        """
        # Convert sequence to complex numbers
        complex_seq = []

        for aa in sequence.upper():
            try:
                complex_seq.append(aa_to_complex[aa])
            except KeyError:
                # Handle unknown amino acids
                complex_seq.append(aa_to_complex['X'])
        # Apply FFT
        complex_seq = np.array(complex_seq)
        
        fft_result = fft(complex_seq, n=2048)
        return fft_result
