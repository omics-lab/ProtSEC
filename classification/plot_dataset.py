import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import sys
from sklearn.manifold import TSNE
import umap


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from embedder import ProteinEmbedder

# Define the directories containing the protein FASTA files
protein_dirs = [
    'ATP synthase',
    'Alpha-amylase',
    'Beta-galactosidase',
    'DNA polymerase I',
    'Protease'
]

# Function to read protein sequences from FASTA files
def read_protein_sequences(directory):
    sequences = []
    directory = os.path.join('dataset', directory)
    # Check if directory exists
    if not os.path.exists(directory):
        print(f"Directory not found: {directory}")
        return sequences
    
    # Get all FASTA files in the directory
    fasta_files = [f for f in os.listdir(directory) if f.endswith(('.fasta', '.fa'))]
    
    for fasta_file in fasta_files:
        file_path = os.path.join(directory, fasta_file)
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                sequences.append(str(record.seq))
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
        
    
    return sequences

def embedding(sequence):
    emb = ProteinEmbedder(n = 768)
    # Generate the embedding for the sequence
    complex_vector = emb.encode(sequence)
    # Convert to a 2D array
    real_vector = np.c_[complex_vector.real, complex_vector.imag].flatten()
    
    return real_vector.tolist()

# Collect all sequences and their corresponding proteins
all_sequences = []
protein_labels = []

for protein_dir in protein_dirs:
    sequences = read_protein_sequences(protein_dir)
    print(f"Found {len(sequences)} sequences in {protein_dir}")
    
    all_sequences.extend(sequences)
    protein_labels.extend([protein_dir] * len(sequences))

# Generate embeddings for all sequences
print("Generating embeddings...")
X = np.array([embedding(seq) for seq in all_sequences])

# # Standardize the features
# scaler = StandardScaler()
# X_scaled = scaler.fit_transform(X)

# Apply PCA to reduce dimensionality to 2
# print("Applying PCA...")
# pca = PCA(n_components=2, random_state=42)
# X_pca = pca.fit_transform(X)

# X_tsne = TSNE(n_components=2, learning_rate='auto', init='random').fit_transform(X)

# Apply UMAP to reduce dimensionality to 2
# print("Applying UMAP...")
reducer = umap.UMAP(
    n_neighbors=100,
    min_dist=0.7,
    n_components=2
)

X_umap = reducer.fit_transform(X)


# Create a 2D scatter plot
plt.figure(figsize=(12, 8))
sns.set_style("whitegrid")

# Create a color palette
palette = sns.color_palette("husl", len(protein_dirs))

# Create a dictionary to map protein names to their respective colors
color_dict = dict(zip(protein_dirs, palette))

# Create a scatter plot with different colors for each protein type
for protein in protein_dirs:
    idx = [i for i, label in enumerate(protein_labels) if label == protein]
    # print(idx)

    if idx:
        plt.scatter(X_umap[idx, 0], X_umap[idx, 1], label=protein, alpha=0.7, s=5)

plt.title('Scatter Plot of Protein Sequences', fontsize=16)
plt.xlabel(f'Component 1', fontsize=12)
plt.ylabel(f'Component 2', fontsize=12)
plt.legend(fontsize=10)
plt.tight_layout()

plt.savefig('protein_clusters_umap.png', dpi=300, bbox_inches='tight')
plt.show()

