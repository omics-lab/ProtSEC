#!/bin/bash

# k-NN Classifier Usage Examples for EC-CARE Task 1
# Make sure to have the required dependencies installed:
# pip install scikit-learn pandas numpy matplotlib seaborn biopython

echo "============================================"
echo "k-NN Classifier for Enzyme Classification"
echo "============================================"

# Make the Python script executable
chmod +x knn.py

# Example 1: EC1 level classification
echo "Running EC1 classification..."
python knn.py \
    --train_matrix ./fastas/all_test_seqs_score_matrix.csv \
    --train_annotations protein_train.csv \
    --test_matrix ./fastas/30_protein_test_score_matrix.csv \
    --test_annotations 30_protein_test.csv \
    --ec_level EC1 \
    --output_dir results_ec1

echo "EC1 classification completed!"

# Example 2: EC2 level classification  
echo "Running EC2 classification..."
python knn.py \
    --train_matrix ./fastas/train_swissprot_1000_ProtSEC_matrix.csv \
    --train_annotations protein_train.csv \
    --test_matrix ./fastas/30_protein_test_score_matrix.csv \
    --test_annotations 30_protein_test.csv \
    --ec_level EC2 \
    --output_dir results_ec2

echo "EC2 classification completed!"

# Example 3: EC3 level classification
echo "Running EC3 classification..."
python knn.py \
    --train_matrix ./fastas/all_test_seqs_score_matrix.csv \
    --train_annotations protein_train.csv \
    --test_matrix ./fastas/30_protein_test_score_matrix.csv \
    --test_annotations 30_protein_test.csv \
    --ec_level EC3 \
    --output_dir results_ec3

echo "EC3 classification completed!"

# Example 4: Full EC number classification
echo "Running full EC number classification..."
python knn.py \
    --train_matrix ./fastas/all_test_seqs_score_matrix.csv \
    --train_annotations protein_train.csv \
    --test_matrix ./fastas/30_protein_test_score_matrix.csv \
    --test_annotations 30_protein_test.csv \
    --ec_level "EC All" \
    --output_dir results_ec_all

echo "Full EC number classification completed!"

# Example 5: Custom k values
echo "Running classification with custom k values..."
python knn.py \
    --train_matrix ./fastas/all_test_seqs_score_matrix.csv \
    --train_annotations protein_train.csv \
    --test_matrix ./fastas/30_protein_test_score_matrix.csv \
    --test_annotations 30_protein_test.csv \
    --ec_level EC1 \
    --k_values 1 3 5 7 9 11 15 \
    --output_dir results_custom_k

echo "Custom k values classification completed!"

echo "============================================"
echo "All classifications completed!"
echo "Check the results_* directories for outputs"
echo "============================================"
