"""
A module for encoding protein sequences using complex number representations.
Main encoding process:
1. Uses BLOSUM62 substitution matrix for amino acid relationships
2. Projects amino acids onto a 2D complex plane using multiple methods (MDS, t-SNE, UMAP)
UMAP(metric='precomputed').fit_transform(<distance matrix>)
"""

import numpy as np
import umap
from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio.Align import substitution_matrices
import warnings
warnings.filterwarnings('ignore')

class AminoAcidEncoder:
    """
    A class to encode amino acids into complex number representations.
    The encoding is based on the BLOSUM62 substitution matrix and uses
    dimensionality reduction techniques to project the amino acids into a
    2D complex plane.
    """

    def __init__(self):
        """
        Initializes the encoder with amino acid properties and dissimilarity matrix.
        """
        self.amino_acids = list("ABCDEFGHIKLMNPQRSTVWYZ")  # Extended amino acid alphabet
        # Alternative standard 20 amino acids: "ARNDCQEGHILKMFPSTWYV"
        self.aa_to_complex = self._create_complex_embeddings()

    def _create_dissimilarity_matrix(self):
        """
        Creates a dissimilarity matrix from the BLOSUM62 substitution matrix.
        
        Returns:
            numpy.ndarray: Dissimilarity matrix
        """
        # Load BLOSUM62 matrix for amino acid similarity scores
        blosum62 = substitution_matrices.load("BLOSUM62")
        
        # Initialize similarity matrix for all amino acid pairs
        n_aa = len(self.amino_acids)
        sim_matrix = np.zeros((n_aa, n_aa))
        
        # Populate similarity matrix with BLOSUM62 scores
        for i, aa1 in enumerate(self.amino_acids):
            for j, aa2 in enumerate(self.amino_acids):
                try:
                    sim_matrix[i, j] = blosum62[(aa1, aa2)]
                except KeyError:
                    # Handle cases where certain amino acids might not be in BLOSUM62
                    sim_matrix[i, j] = -4  # Default value for unknown pairs
        
        # Convert similarity scores to dissimilarity measures
        max_sim = np.max(sim_matrix)
        diss_matrix = max_sim - sim_matrix
        
        return diss_matrix

    def _create_complex_embeddings(self):
        """
        Creates complex number embeddings for amino acids.
        
        Process:
        1. Extracts similarity scores from BLOSUM62 matrix
        2. Converts similarities to dissimilarities
        3. Uses selected method to project amino acids onto 2D space
        4. Maps coordinates to complex numbers
        
        Returns:
            dict: Mapping of amino acids to their complex number representations
        """
        # Get dissimilarity matrix
        diss_matrix = self._create_dissimilarity_matrix()

        
        reducer = umap.UMAP(
            n_components=2,
            metric='precomputed',
            random_state=42
        )
        coords = reducer.fit_transform(diss_matrix)
        

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

    def visualize_embeddings(self, title_suffix=""):
        """
        Creates a scatter plot showing amino acid positions in the complex plane.
        
        Args:
            title_suffix (str): Additional text to append to the plot title
        """
        plt.figure(figsize=(10, 8))
        for aa, c in self.aa_to_complex.items():
            plt.scatter(c.real, c.imag)
            plt.text(c.real + 0.01, c.imag + 0.01, aa, fontsize=12)

        plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
        plt.grid(True, alpha=0.3)
        title = f'Amino Acid Embeddings in Complex Plane from BLOSUM62 using UMAP'
        if title_suffix:
            title += f" - {title_suffix}"
        plt.title(title)
        plt.xlabel('Real Part')
        plt.ylabel('Imaginary Part')
        plt.show()
        
    
# Main execution block for testing
if __name__ == "__main__":
    encoder = AminoAcidEncoder()
    # encoder.visualize_embeddings()
    print(encoder.get_amino_acid_to_complex())
