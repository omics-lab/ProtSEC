# -*- coding: utf-8 -*-
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from Bio.Align import substitution_matrices
from sklearn.manifold import TSNE
import umap

import warnings
warnings.filterwarnings('ignore')






class AminoAcidEncoder:
    """
    Class to encode amino acids into complex number representations using MDS.
    """
    
    def __init__(self, dim_reduct, distance_method='SMS'):
        """
        Initialize the encoder with a specified distance method.
        
        Args:
            distance_method (str): The distance method to use. 
                Options: 'SMS' (original) or 'heur' (new)
        """
        self.amino_acids = list("ABCDEFGHIKLMNPQRSTVWYZ")  # Extended amino acid alphabet
        # Alternative standard 20 amino acids: "ARNDCQEGHILKMFPSTWYV"
        self.dim_reduct = dim_reduct
        self.distance_method = distance_method
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
        
        # Convert similarity scores to dissimilarity measures based on the selected method
        if self.distance_method == 'SMS':
            # Subtract from Max Similarity
            max_sim = np.max(sim_matrix)
            diss_matrix = max_sim - sim_matrix
        elif self.distance_method == 'ASMP':
            # Average of Self Minus Pairwise
            diss_matrix = np.zeros_like(sim_matrix)
            for i in range(n_aa):
                for j in range(n_aa):
                    diss_matrix[i, j] = (sim_matrix[i, i] + sim_matrix[j, j]) / 2 - sim_matrix[i, j] 
        elif self.distance_method == 'SNN':
            # Shift and Normalize (SNN)
            mn = np.min(sim_matrix)
            mx = -mn + np.max(sim_matrix)
            diss_matrix = np.zeros_like(sim_matrix)
            for i in range(n_aa):
                for j in range(n_aa):
                    sim_prime = -mn + sim_matrix[i, j]
                    diss_matrix[i, j] = 1.0 - (sim_prime / mx)
        
        else:
            raise ValueError(f"Unknown distance method: {self.distance_method}")
        # print(diss_matrix)
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

        if self.dim_reduct == 't-SNE':
            # Use t-SNE for dimensionality reduction
            reducer = TSNE(n_components=2,
                           metric='precomputed',
                           init='random',
                           perplexity=8,
                           random_state=31)
        elif self.dim_reduct == 'UMAP': 
            # Use UMAP for dimensionality reduction
            reducer = umap.UMAP(n_components=2, metric='precomputed', random_state=31)
        elif self.dim_reduct == 'MDS':
            # Use MDS for dimensionality reduction  
            reducer = MDS(n_components=2, dissimilarity='precomputed', random_state=41)
        # Note: The 'precomputed' metric is used because we are passing a dissimilarity matrix
        # Fit and transform the dissimilarity matrix to get 2D coordinates
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
        
        title = f'Amino Acid Embeddings in Complex Plane from BLOSUM62 ({self.dim_reduct} + {self.distance_method})'
        
        plt.title(title)
        plt.xlabel('Real Part')
        plt.ylabel('Imaginary Part')
        # plt.show(block=True)
        plt.savefig('temp.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    
# Main execution block for testing
if __name__ == "__main__":
    
    
    # Create encoder with new distance method
    method, fun = 'UMAP', 'SNN'
    encoder = AminoAcidEncoder(method, distance_method=fun)
    encoder.visualize_embeddings()
    mp = encoder.get_amino_acid_to_complex()
    print(mp)
    