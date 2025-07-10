import numpy as np
from embedder import ProteinEmbedder
import xgboost as xgb
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
import pandas as pd
from Bio import SeqIO
import sys

def cut_protein_sequence_center(protein_sequence, n):
    # If sequence is already shorter than or equal to n, return as is
    if len(protein_sequence) <= n:
        return protein_sequence
    
    # Calculate how many residues to keep from each end
    total_to_keep = n
    
    # Split the kept residues equally between N-terminal and C-terminal
    # If odd number to keep, keep one more from the N-terminal
    n_terminal_keep = (total_to_keep) // 2
    c_terminal_keep = total_to_keep - n_terminal_keep
    
    # Extract N-terminal and C-terminal portions
    n_terminal = protein_sequence[:n_terminal_keep]
    c_terminal = protein_sequence[-c_terminal_keep:] if c_terminal_keep > 0 else ""
    
    # Combine N-terminal and C-terminal portions
    return n_terminal + c_terminal

def plot_confusion_matrix(y_true, y_pred, classes, title='Confusion Matrix'):
    """Plot confusion matrix with proper formatting"""
    cm = confusion_matrix(y_true, y_pred)
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=classes, yticklabels=classes,
                cbar_kws={'label': 'Count'})
    plt.title('MultiLoc 4 Class Confusion Matrix', fontsize=16, fontweight='bold')
    plt.xlabel('Predicted Location', fontsize=12)
    plt.ylabel('True Location', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()

def calculate_metrics(y_true, y_pred, classes):
    """Calculate and display comprehensive classification metrics"""
    # Basic metrics
    accuracy = accuracy_score(y_true, y_pred)
    precision_macro = precision_score(y_true, y_pred, average='macro')
    recall_macro = recall_score(y_true, y_pred, average='macro')
    f1_macro = f1_score(y_true, y_pred, average='macro')
    
    precision_micro = precision_score(y_true, y_pred, average='micro')
    recall_micro = recall_score(y_true, y_pred, average='micro')
    f1_micro = f1_score(y_true, y_pred, average='micro')
    
    precision_weighted = precision_score(y_true, y_pred, average='weighted')
    recall_weighted = recall_score(y_true, y_pred, average='weighted')
    f1_weighted = f1_score(y_true, y_pred, average='weighted')
    
    # Print summary metrics
    print("\n" + "="*60)
    print("CLASSIFICATION METRICS SUMMARY")
    print("="*60)
    print(f"Accuracy: {accuracy:.4f}")
    print(f"\nMacro-averaged metrics:")
    print(f"  Precision: {precision_macro:.4f}")
    print(f"  Recall:    {recall_macro:.4f}")
    print(f"  F1-score:  {f1_macro:.4f}")
    print(f"\nMicro-averaged metrics:")
    print(f"  Precision: {precision_micro:.4f}")
    print(f"  Recall:    {recall_micro:.4f}")
    print(f"  F1-score:  {f1_micro:.4f}")
    print(f"\nWeighted-averaged metrics:")
    print(f"  Precision: {precision_weighted:.4f}")
    print(f"  Recall:    {recall_weighted:.4f}")
    print(f"  F1-score:  {f1_weighted:.4f}")
    
    # Detailed classification report
    print("\n" + "="*60)
    print("DETAILED CLASSIFICATION REPORT")
    print("="*60)
    print(classification_report(y_true, y_pred, target_names=classes, digits=4))
    
    return {
        'accuracy': accuracy,
        'precision_macro': precision_macro,
        'recall_macro': recall_macro,
        'f1_macro': f1_macro,
        'precision_micro': precision_micro,
        'recall_micro': recall_micro,
        'f1_micro': f1_micro,
        'precision_weighted': precision_weighted,
        'recall_weighted': recall_weighted,
        'f1_weighted': f1_weighted
    }

if __name__ == "__main__":
    # Load data
    clf_data_path = 'classification/Multiloc/merged_multiloc.fasta'
    train_data, test_data = [], []
    selected_labels = ['Cell Membrane', 'Cytoplasm', 'Extracellular', 'Nucleus']

    mp = {
     'Cell.membrane-U': 'Cell Membrane',
     'Cytoplasm-U': 'Cytoplasm',
     'Endoplasmic.reticulum-U': 'Endoplasmic Reticulum',
     'Extracellular-U': 'Extracellular',
     'Golgi.apparatus-U': 'Golgi Apparatus',
     'Lysosome/Vacuole-U': 'Lysosome/Vacuole',
     'Mitochondrion-U': 'Mitochondrion',
     'Nucleus-U': 'Nucleus',
     'Peroxisome-U': 'Peroxisome',
     'Plastid-U': 'Plastid'
    }

    for record in SeqIO.parse(clf_data_path, 'fasta'):
        header = str(record.description)
        label = header.split(' ')[1]
        label = mp[label]
        if label not in selected_labels:
            continue
        if 'test' in header:
            test_data.append((str(record.seq), label))
        else:
            train_data.append((str(record.seq), label))
    
    print(f"Loaded {len(train_data)} training samples and {len(test_data)} test samples.")
    
    dim = 768  # Dimension for embedding
    embedder = ProteinEmbedder(dim_reduct='MDS', dist_func='SMS', n=dim)
    
    # Encode protein sequences
    print("Encoding protein sequences...")
    embedded_training_data = []
    for seq, label in train_data:
        v = embedder.encode(cut_protein_sequence_center(seq, dim))
        embedded_training_data.append((v, label))
    
    embedded_test_data = []
    for seq, label in test_data:
        v = embedder.encode(cut_protein_sequence_center(seq, dim))
        embedded_test_data.append((v, label))
    
    np.random.shuffle(embedded_training_data)
    np.random.shuffle(embedded_test_data)

    X_train, y_train = map(list, zip(*embedded_training_data))
    X_test, y_test = map(list, zip(*embedded_test_data))

    X_train = np.array(X_train)
    X_test = np.array(X_test)

    # Convert complex embeddings to real features
    X_train = np.concatenate([X_train.real, X_train.imag], axis=1)
    X_test = np.concatenate([X_test.real, X_test.imag], axis=1)

    # Encode labels
    label_encoder = LabelEncoder()
    y_train_encoded = label_encoder.fit_transform(y_train)
    y_test_encoded = label_encoder.transform(y_test)

    print(f"Number of classes: {len(label_encoder.classes_)}")
    print(f"Class labels: {label_encoder.classes_}")
    print(f"Training set shape: {X_train.shape}")
    print(f"Test set shape: {X_test.shape}")

    # Train XGBoost classifier
    print("\n" + "="*60)
    print("TRAINING XGBOOST CLASSIFIER")
    print("="*60)
    
    # XGBoost parameters
    xgb_params = {
        'objective': 'multi:softprob',
        'num_class': len(label_encoder.classes_),
        'max_depth': 5,
        'learning_rate': 0.1,
        'n_estimators': 1000,
        'subsample': 0.8,
        'colsample_bytree': 0.8,
        'random_state': 37,
        'eval_metric': 'mlogloss',
        'verbosity': 1,
    }
    
    # Create and train XGBoost classifier
    xgb_classifier = xgb.XGBClassifier(**xgb_params)
    
    print("Training XGBoost classifier...")
    xgb_classifier.fit(X_train, y_train_encoded)
    
    # Make predictions
    print("Making predictions...")
    y_pred_encoded = xgb_classifier.predict(X_test)
    y_pred_proba = xgb_classifier.predict_proba(X_test)
    
    # Calculate and display metrics
    metrics = calculate_metrics(y_test_encoded, y_pred_encoded, label_encoder.classes_)
    
    # Plot confusion matrix
    print("\nGenerating confusion matrix...")
    plot_confusion_matrix(y_test_encoded, y_pred_encoded, 
                         label_encoder.classes_, 
                         title='XGBoost Classifier - Confusion Matrix')
    
    """
    ['Cell.membrane-U' 'Cytoplasm-U' 'Endoplasmic.reticulum-U'
 'Extracellular-U' 'Golgi.apparatus-U' 'Lysosome/Vacuole-U'
 'Mitochondrion-U' 'Nucleus-U' 'Peroxisome-U' 'Plastid-U']

 mp = {
     'Cell.membrane-U': 'Cell Membrane',
     'Cytoplasm-U': 'Cytoplasm',
     'Endoplasmic.reticulum-U': 'Endoplasmic Reticulum',
     'Extracellular-U': 'Extracellular',
     'Golgi.apparatus-U': 'Golgi Apparatus',
     'Lysosome/Vacuole-U': 'Lysosome/Vacuole',
     'Mitochondrion-U': 'Mitochondrion',
     'Nucleus-U': 'Nucleus',
     'Peroxisome-U': 'Peroxisome',
     'Plastid-U': 'Plastid'
 }
    """