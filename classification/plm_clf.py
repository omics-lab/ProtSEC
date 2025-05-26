import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
import seaborn as sns
import sys
from sklearn.manifold import TSNE
from sklearn.metrics import confusion_matrix, classification_report, accuracy_score
from sklearn.neural_network import MLPClassifier
import xgboost as xgb
import pandas as pd

# Define the directories containing the protein embedding files
protein_dirs = [
    'ATP synthase',
    'Alpha-amylase',
    'Beta-galactosidase', 
    'DNA polymerase I',
    'Protease'
]

# Base directory for ESM2 embeddings
esm2_base_dir = 'clf_vector/prot_bert'

# Function to load embeddings from .npy files
def load_embeddings(directory, file_name):
    embeddings = []
    full_path = os.path.join(esm2_base_dir, directory, file_name)
    
    # Check if file exists
    if not os.path.exists(full_path):
        print(f"File not found: {full_path}")
        return np.array([])
    
    try:
        embeddings = np.load(full_path)
        print(f"Loaded {len(embeddings)} embeddings from {full_path}")
        print(f"Embedding shape: {embeddings.shape}")
    except Exception as e:
        print(f"Error loading {full_path}: {e}")
        return np.array([])
    
    return embeddings

# Collect training embeddings and their corresponding protein labels
X_train_list = []
y_train_list = []

print("Loading training embeddings...")
for protein_dir in protein_dirs:
    embeddings = load_embeddings(protein_dir, 'train_data.npy')
    if len(embeddings) > 0:
        X_train_list.append(embeddings)
        y_train_list.extend([protein_dir] * len(embeddings))

# Collect test embeddings and their corresponding protein labels
X_test_list = []
y_test_list = []

print("Loading test embeddings...")
for protein_dir in protein_dirs:
    embeddings = load_embeddings(protein_dir, 'test_data.npy')
    if len(embeddings) > 0:
        X_test_list.append(embeddings)
        y_test_list.extend([protein_dir] * len(embeddings))

# Concatenate all embeddings
if X_train_list:
    X_train = np.vstack(X_train_list)
    y_train = np.array(y_train_list)
else:
    print("No training data loaded!")
    sys.exit(1)

if X_test_list:
    X_test = np.vstack(X_test_list)
    y_test = np.array(y_test_list)
else:
    print("No test data loaded!")
    sys.exit(1)

print(f"Training embeddings shape: {X_train.shape}")
print(f"Test embeddings shape: {X_test.shape}")

# Optional: Standardize the features (uncomment if needed)
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
