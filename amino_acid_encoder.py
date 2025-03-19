"""
A module for encoding protein sequences using complex number representations.
Main encoding process:
1. Uses BLOSUM62 substitution matrix for amino acid relationships
2. Projects amino acids onto a 2D complex plane
"""

import numpy as np
from sklearn.manifold import MDS
from scipy.fftpack import fft
from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio.Align import substitution_matrices
import warnings
warnings.filterwarnings('ignore')

class AminoAcidEncoder:
    """
    A class for encoding protein sequences into fixed-length numerical vectors.
    
    The encoding process uses:
    - Complex number embeddings derived from BLOSUM62 matrix
    - Multidimensional scaling (MDS) for dimensionality reduction
    - Fast Fourier Transform for sequence processing
    """

    def __init__(self):
        """
        Initializes the encoder with default amino acid alphabet.
        Extended alphabet includes B, X, Z for handling ambiguous amino acids.
        """
        self.amino_acids = list("ABCDEFGHIKLMNPQRSTVWXYZ")  # Extended amino acid alphabet
        # Alternative standard 20 amino acids: "ARNDCQEGHILKMFPSTWYV"
        self.aa_to_complex = self._create_complex_embeddings()

    def _create_complex_embeddings(self):
        """
        Creates complex number embeddings for amino acids.
        
        Process:
        1. Extracts similarity scores from BLOSUM62 matrix
        2. Converts similarities to dissimilarities
        3. Uses MDS to project amino acids onto 2D space
        4. Maps coordinates to complex numbers
        
        Returns:
            dict: Mapping of amino acids to their complex number representations
        """
        # Load BLOSUM62 matrix for amino acid similarity scores
        blosum62 = substitution_matrices.load("BLOSUM62")
        
        # Initialize similarity matrix for all amino acid pairs
        n_aa = len(self.amino_acids)
        sim_matrix = np.zeros((n_aa, n_aa))
        
        # Populate similarity matrix with BLOSUM62 scores
        for i, aa1 in enumerate(self.amino_acids):
            for j, aa2 in enumerate(self.amino_acids):
                sim_matrix[i, j] = blosum62[(aa1, aa2)]
        
        # Convert similarity scores to dissimilarity measures
        max_sim = np.max(sim_matrix)
        diss_matrix = max_sim - sim_matrix

        # Project amino acids onto 2D space using MDS
        mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
        coords = mds.fit_transform(diss_matrix)

        # Scale coordinates to [0,1] range
        coords = (coords - np.min(coords)) / (np.max(coords) - np.min(coords))

        # Create mapping from amino acids to complex numbers
        aa_to_complex = {}
        for i, aa in enumerate(self.amino_acids):
            aa_to_complex[aa] = complex(coords[i, 0], coords[i, 1])

        return aa_to_complex

    def get_amino_acid_to_complex(self):
        """Returns the dictionary mapping amino acids to their complex representations."""
        return self.aa_to_complex

    def visualize_embeddings(self):
        """
        Creates a scatter plot showing amino acid positions in the complex plane.
        
        Plot features:
        - Each amino acid is represented as a point
        - Axes show real and imaginary components
        - Grid lines and labels for better visualization
        """
        plt.figure(figsize=(10, 8))
        for aa, c in self.aa_to_complex.items():
            plt.scatter(c.real, c.imag)
            plt.text(c.real + 0.01, c.imag + 0.01, aa, fontsize=12)

        plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
        plt.grid(True, alpha=0.3)
        plt.title('Amino Acid Embeddings in Complex Plane from BLOSUM62')
        plt.xlabel('Real Part')
        plt.ylabel('Imaginary Part')
        plt.show()

# Main execution block for testing
if __name__ == "__main__":
    amino_acid_encoder = AminoAcidEncoder()
    # Uncomment to visualize the embeddings
    # amino_acid_encoder.visualize_embeddings()
    amino_acid_to_complex = amino_acid_encoder.get_amino_acid_to_complex()
    print(amino_acid_to_complex)
