import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
import seaborn as sns
import sys
from sklearn.manifold import TSNE
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, classification_report, accuracy_score
from sklearn.neural_network import MLPClassifier
import xgboost as xgb
import pandas as pd

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
    file_paths = []
    directory = os.path.join('dataset', directory)
    
    # Check if directory exists
    if not os.path.exists(directory):
        print(f"Directory not found: {directory}")
        return sequences, file_paths
    
    # Get all FASTA files in the directory
    fasta_files = [f for f in os.listdir(directory) if f.endswith(('.fasta', '.fa'))]
    
    for fasta_file in fasta_files:
        file_path = os.path.join(directory, fasta_file)
        file_sequences = []
        
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                file_sequences.append(str(record.seq))
            
            # Add sequences and their source file
            sequences.extend(file_sequences)
            file_paths.extend([fasta_file] * len(file_sequences))
            
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
    
    return sequences, file_paths

def embedding(sequence):
    emb = ProteinEmbedder(n=768)
    # Generate the embedding for the sequence
    complex_vector = emb.encode(sequence)
    # Convert to a 2D array
    real_vector = np.c_[complex_vector.real, complex_vector.imag].flatten()
    
    return real_vector.tolist()

# Collect all sequences, their corresponding proteins, and source files
all_sequences = []
protein_labels = []
file_sources = []

for protein_dir in protein_dirs:
    sequences, file_paths = read_protein_sequences(protein_dir)
    print(f"Found {len(sequences)} sequences in {protein_dir}")
    
    all_sequences.extend(sequences)
    protein_labels.extend([protein_dir] * len(sequences))
    file_sources.extend(file_paths)

# Generate embeddings for all sequences
print("Generating embeddings...")
X = np.array([embedding(seq) for seq in all_sequences])
y = np.array(protein_labels)

# Create a dataframe to keep track of sequences, labels, and file sources
data_df = pd.DataFrame({
    'sequence': all_sequences,
    'label': protein_labels,
    'file_source': file_sources
})

# Create train/test splits with 5% from each fasta file selected for testing
X_train_indices = []
X_test_indices = []

# Group by protein type and file source
for protein in protein_dirs:
    protein_files = data_df[data_df['label'] == protein]['file_source'].unique()
    
    for file in protein_files:
        file_indices = data_df[(data_df['label'] == protein) & 
                               (data_df['file_source'] == file)].index.tolist()
        
        if len(file_indices) > 0:
            # Split indices into train (95%) and test (5%) sets
            test_size = max(1, int(len(file_indices) * 0.05))  # Ensure at least one test sample
            test_indices = np.random.choice(file_indices, size=test_size, replace=False)
            train_indices = list(set(file_indices) - set(test_indices))
            
            X_train_indices.extend(train_indices)
            X_test_indices.extend(test_indices)

# Create the actual train and test sets
X_train = X[X_train_indices]
X_test = X[X_test_indices]
y_train = y[X_train_indices]
y_test = y[X_test_indices]

# Standardize the features
scaler = StandardScaler()
# X_train_scaled = scaler.fit_transform(X_train)
# X_test_scaled = scaler.transform(X_test)

# Encode string labels to integers using LabelEncoder
label_encoder = LabelEncoder()
y_train_encoded = label_encoder.fit_transform(y_train)
y_test_encoded = label_encoder.transform(y_test)

# Store the original class names for later use
class_names = label_encoder.classes_
print(f"Encoded classes: {list(zip(class_names, range(len(class_names))))}")

print(f"Training set size: {len(X_train)}")
print(f"Test set size: {len(X_test)}")

print("\nTraining XGBoost classifier...")
# Train XGBoost classifier
xgb_clf = xgb.XGBClassifier(
    objective='multi:softprob',
    n_estimators=100,
    max_depth=5,
    learning_rate=0.1,
    random_state=42
)
xgb_clf.fit(X_train, y_train_encoded)

# Evaluate XGBoost
xgb_preds_encoded = xgb_clf.predict(X_test)
xgb_preds = label_encoder.inverse_transform(xgb_preds_encoded)
xgb_accuracy = accuracy_score(y_test, xgb_preds)
xgb_confusion = confusion_matrix(y_test, xgb_preds)
xgb_report = classification_report(y_test, xgb_preds)

print(f"XGBoost Accuracy: {xgb_accuracy:.4f}")
print("\nXGBoost Confusion Matrix:")
print(xgb_confusion)
print("\nXGBoost Classification Report:")
print(xgb_report)

print("\nTraining MLP classifier...")
# Train MLP classifier
mlp_clf = MLPClassifier(
    hidden_layer_sizes=(100, 50),
    activation='relu',
    solver='adam',
    alpha=0.0001,
    max_iter=300,
    random_state=42,
    early_stopping=True
)
mlp_clf.fit(X_train, y_train_encoded)

# Evaluate MLP
mlp_preds_encoded = mlp_clf.predict(X_test)
mlp_preds = label_encoder.inverse_transform(mlp_preds_encoded)
mlp_accuracy = accuracy_score(y_test, mlp_preds)
mlp_confusion = confusion_matrix(y_test, mlp_preds)
mlp_report = classification_report(y_test, mlp_preds)

print(f"MLP Accuracy: {mlp_accuracy:.4f}")
print("\nMLP Confusion Matrix:")
print(mlp_confusion)
print("\nMLP Classification Report:")
print(mlp_report)

# Plot confusion matrices
def plot_confusion_matrix(cm, classes, title, normalize=False, cmap=plt.cm.Blues):
    """
    Plot confusion matrix.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        fmt = '.2f'
    else:
        fmt = 'd'

    plt.figure(figsize=(10, 8))
    sns.heatmap(cm, annot=True, fmt=fmt, cmap=cmap, 
                xticklabels=classes, yticklabels=classes)
    plt.title(title, fontsize=14)
    plt.ylabel('True label', fontsize=12)
    plt.xlabel('Predicted label', fontsize=12)
    plt.tight_layout()
    plt.savefig(f'{title.replace(" ", "_").lower()}.png', dpi=300, bbox_inches='tight')
    plt.show()

# Plot XGBoost confusion matrix
plot_confusion_matrix(
    xgb_confusion, 
    classes=protein_dirs,
    title='XGBoost Confusion Matrix'
)

# Plot MLP confusion matrix
plot_confusion_matrix(
    mlp_confusion, 
    classes=protein_dirs,
    title='MLP Confusion Matrix'
)

# Also visualize the test data with TSNE
X_test_tsne = TSNE(n_components=2, learning_rate='auto', init='random', random_state=42).fit_transform(X_test)

# Create a 2D scatter plot of test data
plt.figure(figsize=(12, 8))
sns.set_style("whitegrid")

# Create a color palette
palette = sns.color_palette("husl", len(protein_dirs))

# Create a dictionary to map protein names to their respective colors
color_dict = dict(zip(protein_dirs, palette))

# Create a scatter plot with different colors for each protein type
for protein in protein_dirs:
    idx = [i for i, label in enumerate(y_test) if label == protein]
    if idx:
        plt.scatter(X_test_tsne[idx, 0], X_test_tsne[idx, 1], label=protein, alpha=0.7, s=30)

plt.title('TSNE of Test Protein Sequences', fontsize=16)
plt.xlabel('Component 1', fontsize=12)
plt.ylabel('Component 2', fontsize=12)
plt.legend(fontsize=10)
plt.tight_layout()

plt.savefig('test_protein_clusters_tsne.png', dpi=300, bbox_inches='tight')
plt.show()
