import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
import seaborn as sns
import sys
from sklearn.manifold import TSNE
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
def read_protein_sequences(directory, file_name):
    sequences = []
    directory = os.path.join('dataset', directory)
    
    # Check if directory exists
    if not os.path.exists(directory):
        print(f"Directory not found: {directory}")
        return sequences
    
    # Check if specific file exists
    file_path = os.path.join(directory, file_name)
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return sequences
    
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            sequences.append(str(record.seq))
        print(f"Read {len(sequences)} sequences from {file_path}")
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    
    return sequences

def embedding(sequence):
    emb = ProteinEmbedder(n=768)
    # Generate the embedding for the sequence
    complex_vector = emb.encode(sequence)
    # Convert to a 2D array
    real_vector = np.c_[complex_vector.real, complex_vector.imag].flatten()
    
    return real_vector.tolist()

# Collect training sequences and their corresponding protein labels
train_sequences = []
train_labels = []

print("Reading training data...")
for protein_dir in protein_dirs:
    sequences = read_protein_sequences(protein_dir, 'train_data.fasta')
    train_sequences.extend(sequences)
    train_labels.extend([protein_dir] * len(sequences))

# Collect test sequences and their corresponding protein labels
test_sequences = []
test_labels = []

print("Reading testing data...")
for protein_dir in protein_dirs:
    sequences = read_protein_sequences(protein_dir, 'test_data.fasta')
    test_sequences.extend(sequences)
    test_labels.extend([protein_dir] * len(sequences))

# Generate embeddings for training sequences
print("Generating embeddings for training data...")
X_train = np.array([embedding(seq) for seq in train_sequences])
y_train = np.array(train_labels)

# Generate embeddings for test sequences
print("Generating embeddings for testing data...")
X_test = np.array([embedding(seq) for seq in test_sequences])
y_test = np.array(test_labels)

# Standardize the features
# scaler = StandardScaler()
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
    title='Confusion Matrix (ProSEC)'
)

