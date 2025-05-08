from Bio import SeqIO
from merged_class import sub_to_super
import numpy as np
import sys
import os
import xgboost as xgb
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix
from sklearn.decomposition import PCA
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns


import random
random.seed(31)

# Add the parent directory (project/) to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from embedder import ProteinEmbedder

def normalize(x: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(x)
    return x / (norm if norm > 0 else 1)



def get_embedding_cls_pairs(data_path):
    """
    Get the embedding and class pairs from the data path.
    """
    emb = ProteinEmbedder(n = 1024)

    # unique_cls_to_int = {'Cytoplasm': 0, 'Cell.membrane': 1, 'Nucleus': 2, 'Lysosome/Vacuole': 3, 'Peroxisome': 4, 'Extracellular': 5, 'Endoplasmic.reticulum': 6, 'Plastid': 7, 'Golgi.apparatus': 8, 'Mitochondrion': 9}
    unique_cls_to_int = {'Cytoplasm': 0, 'Cell.membrane': 1, 'Nucleus': 2, 'Extracellular': 3, 'Mitochondrion': 4}
    int_to_cls = {v: k for k, v in unique_cls_to_int.items()}  # For later use with confusion matrix


    list_of_emb_class_pairs = []
    # Read the FASTA file

    for record in SeqIO.parse(data_path, "fasta"):
        sub_class = record.description.split(" ")[1]
        key = sub_to_super[sub_class]
        if key not in unique_cls_to_int:
            continue
        super_class = unique_cls_to_int[key]

        seq = str(record.seq)

        complex_vector = emb.encode(seq)
        real_vector =  np.c_[complex_vector.real, complex_vector.imag].flatten()
        # Normalize the vector
        
        list_of_emb_class_pairs.append((real_vector.tolist(), super_class))
    
    return list_of_emb_class_pairs, int_to_cls



def main():
    train_data_path = "../data/deeplock/train_deeploc_data.fasta"
    test_data_path = "../data/deeplock/test_deeploc_data.fasta"

    train_emb_class_pairs, int_to_cls = get_embedding_cls_pairs(train_data_path)
    test_emb_class_pairs, _ = get_embedding_cls_pairs(test_data_path)

    X_train =  np.array([pair[0] for pair in train_emb_class_pairs])
    y_train = np.array([pair[1] for pair in train_emb_class_pairs])

    # mp = defaultdict(int)
    # for i in y_train:
    #     mp[i] += 1
    # print(mp)

    X_test =  np.array([pair[0] for pair in test_emb_class_pairs])
    y_test = np.array([pair[1] for pair in test_emb_class_pairs])


    pca = PCA(n_components=64, random_state=42)
    X_train_pca = pca.fit_transform(X_train)
    X_test_pca = pca.transform(X_test)

    # Convert to DMatrix for XGBoost
    dtrain = xgb.DMatrix(X_train_pca, label=y_train)
    dval = xgb.DMatrix(X_test_pca, label=y_test)
    # Set parameters for XGBoost
    params = {
        'objective': 'multi:softmax',
        'num_class': 5,
        'max_depth': 5,
        'eta': 0.1,
        'eval_metric': 'mlogloss',
        'subsample': 0.9
    }

    # Train the model
    num_rounds = 100
    evallist = [(dtrain, 'train'), (dval, 'eval')]
    
    model = xgb.train(
        params,
        dtrain,
        num_rounds,
        evallist,
        early_stopping_rounds=20,
        verbose_eval=20
    )

    print("Evaluating model on validation set...")
    y_pred = model.predict(dval)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"Validation accuracy: {accuracy:.4f}")

    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))

    # Generate confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    
    # Create labels for the confusion matrix
    class_names = [int_to_cls[i] for i in range(5)]
    
    # Plot confusion matrix
    plt.figure(figsize=(12, 10))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=class_names, yticklabels=class_names)
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title('Confusion Matrix for Protein Localization')
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    
    # Save confusion matrix plot
    plt.tight_layout()
    plt.savefig('protein_localization_confusion_matrix.png', dpi=300)
    plt.show()

    # Print some additional analysis of the confusion matrix
    print("\nConfusion Matrix Analysis:")
    for i in range(len(class_names)):
        true_positives = cm[i, i]
        false_negatives = cm[i, :].sum() - true_positives
        false_positives = cm[:, i].sum() - true_positives
        
        precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
        recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
        
        print(f"{class_names[i]}:")
        print(f"  - Precision: {precision:.4f}")
        print(f"  - Recall: {recall:.4f}")
        
        # Find most common misclassifications for this class
        if false_negatives > 0:
            wrong_predictions = [(class_names[j], cm[i, j]) for j in range(len(class_names)) if j != i and cm[i, j] > 0]
            if wrong_predictions:
                wrong_predictions.sort(key=lambda x: x[1], reverse=True)
                print(f"  - Most confused with: {wrong_predictions[0][0]} ({wrong_predictions[0][1]} instances)")




if __name__ == "__main__":
    main()
