# k-NN Classifier for Enzyme Class Prediction

This directory contains a comprehensive k-NN (k-Nearest Neighbors) classifier implementation for predicting enzyme classes using ProtSEC-generated distance matrices.

## Features

### 🔬 **Multi-level Classification**
- **EC1**: Main enzyme classes (1-7)
- **EC2**: Subclasses (e.g., 1.1, 1.2)
- **EC3**: Sub-subclasses (e.g., 1.1.1, 1.1.2)
- **EC All**: Complete EC numbers (e.g., 1.1.1.35)

### 🎯 **Advanced k-NN Implementation**
- **Distance Matrix Input**: Uses precomputed ProtSEC distance matrices
- **Similarity Conversion**: Automatically converts distances to similarities (1-distance)
- **Multiple k Values**: Tests different k values to find optimal performance
- **Cross-Validation**: 3-fold CV for robust performance estimation
- **Comprehensive Metrics**: Accuracy, Precision, Recall, F1-score

### 📊 **Visualization & Analysis**
- Performance plots showing accuracy vs k values
- Cross-validation confidence intervals
- Classification reports with per-class metrics
- Prediction confidence scores

### 💾 **Output Files**
- **Predictions**: CSV files with protein IDs and predicted classes
- **Detailed Results**: Text files with comprehensive performance metrics
- **Performance Plots**: High-quality visualizations
- **Model Parameters**: Best k values and hyperparameters

## File Structure

```
📁 analysis/EC_CARE_task1/
├── 📄 knn.py                    # Main k-NN classifier implementation
├── 📄 run_knn_examples.sh       # Usage examples script
├── 📄 README.md                 # This documentation
├── 📄 protein_train.csv         # Training annotations
├── 📄 30_protein_test.csv       # Test annotations
├── 📄 price_protein_test.csv    # Additional test set
└── 📁 fastas/
    ├── 📄 all_test_seqs_score_matrix.csv     # Training distance matrix
    ├── 📄 30_protein_test_score_matrix.csv   # Test distance matrix
    └── 📄 ...                               # Other distance matrices
```

## Installation

### Required Dependencies

```bash
# Essential packages
pip install scikit-learn pandas numpy matplotlib seaborn biopython

# Or install from project requirements
pip install -r ../../requirements.txt
```

### Package Versions (Recommended)
- `scikit-learn >= 1.0.0`
- `pandas >= 1.3.0`
- `numpy >= 1.21.0`
- `matplotlib >= 3.5.0`
- `seaborn >= 0.11.0`
- `biopython >= 1.79`

## Usage

### Basic Usage

```bash
python knn.py \
    --train_matrix ./fastas/all_test_seqs_score_matrix.csv \
    --train_annotations protein_train.csv \
    --test_matrix ./fastas/30_protein_test_score_matrix.csv \
    --test_annotations 30_protein_test.csv \
    --ec_level EC1
```

### Command Line Arguments

| Argument | Required | Description | Example |
|----------|----------|-------------|---------|
| `--train_matrix` | ✅ | Training distance matrix CSV | `./fastas/all_test_seqs_score_matrix.csv` |
| `--train_annotations` | ✅ | Training annotations CSV | `protein_train.csv` |
| `--test_matrix` | ✅ | Test distance matrix CSV | `./fastas/30_protein_test_score_matrix.csv` |
| `--test_annotations` | ✅ | Test annotations CSV | `30_protein_test.csv` |
| `--ec_level` | ❌ | EC level to predict | `EC1` (default), `EC2`, `EC3`, `"EC All"` |
| `--k_values` | ❌ | k values to test | `3 5 7 9 11` (default) |
| `--output_dir` | ❌ | Output directory | `knn_results` (default) |

### Example Commands

#### 1. EC1 Level Classification
```bash
python knn.py \
    --train_matrix ./fastas/all_test_seqs_score_matrix.csv \
    --train_annotations protein_train.csv \
    --test_matrix ./fastas/30_protein_test_score_matrix.csv \
    --test_annotations 30_protein_test.csv \
    --ec_level EC1 \
    --output_dir results_ec1
```

#### 2. Full EC Number Classification
```bash
python knn.py \
    --train_matrix ./fastas/all_test_seqs_score_matrix.csv \
    --train_annotations protein_train.csv \
    --test_matrix ./fastas/30_protein_test_score_matrix.csv \
    --test_annotations 30_protein_test.csv \
    --ec_level "EC All" \
    --output_dir results_ec_all
```

#### 3. Custom k Values
```bash
python knn.py \
    --train_matrix ./fastas/all_test_seqs_score_matrix.csv \
    --train_annotations protein_train.csv \
    --test_matrix ./fastas/30_protein_test_score_matrix.csv \
    --test_annotations 30_protein_test.csv \
    --ec_level EC2 \
    --k_values 1 3 5 7 9 11 15 20 \
    --output_dir results_custom_k
```

### Batch Processing

Run all examples with different EC levels:

```bash
chmod +x run_knn_examples.sh
./run_knn_examples.sh
```

## Input Data Format

### Distance Matrix CSV Format
```csv
ID,0,1,2,3,4,...
0,0.0,0.799,0.865,0.869,...
1,0.799,0.0,0.879,0.884,...
2,0.865,0.879,0.0,0.882,...
...
```

### Annotations CSV Format
```csv
Entry,Entry Name,Sequence,EC number,Length,EC All,clusterRes30,clusterRes50,clusterRes70,clusterRes90,EC3,EC2,EC1
O32178,FADN_BACSU,MHKH...,1.1.1.35,789,1.1.1.35,O32178,O32178,O32178,O32178,1.1.1,1.1,1
P40580,BZRD_YEAST,MGKV...,1.1.1.320,263,1.1.1.320,P40580,P40580,P40580,P40580,1.1.1,1.1,1
...
```

## Output Files

### 1. Predictions CSV
```csv
Protein_ID,Predicted_EC1,Confidence
O32178,1,0.85
P40580,1,0.92
...
```

### 2. Detailed Results Text
```
k-NN Classifier Results for EC1
==================================================

Training data: ./fastas/all_test_seqs_score_matrix.csv
Test data: ./fastas/30_protein_test_score_matrix.csv
EC level: EC1

k=3:
  Accuracy: 0.8500
  Precision: 0.8456
  Recall: 0.8500
  F1-score: 0.8423
  CV Score: 0.8234 ± 0.0234

k=5:
  Accuracy: 0.8700
  Precision: 0.8689
  Recall: 0.8700
  F1-score: 0.8654
  CV Score: 0.8456 ± 0.0198

Best k: 5

Classification Report (Best Model):
              precision    recall  f1-score   support
           1       0.89      0.85      0.87        20
           2       0.83      0.91      0.87        11
           3       0.92      0.85      0.88        13
    accuracy                           0.87        44
   macro avg       0.88      0.87      0.87        44
weighted avg       0.87      0.87      0.87        44
```

### 3. Performance Plots
- **Accuracy vs k**: Shows test and cross-validation accuracy
- **F1-score vs k**: Performance metric across different k values
- **Best Model Summary**: Key performance metrics

## Algorithm Details

### 1. **Data Preprocessing**
```python
# Convert distance to similarity
similarity = 1 - distance_matrix

# Ensure valid similarities
similarity = np.clip(similarity, 0, 1)
np.fill_diagonal(similarity, 1.0)
```

### 2. **k-NN Classification**
```python
# Use precomputed distance matrix
knn = KNeighborsClassifier(n_neighbors=k, metric='precomputed')

# Convert similarity back to distance for sklearn
distance_matrix = 1 - similarity_matrix
knn.fit(distance_matrix, labels)
```

### 3. **Performance Evaluation**
- **Accuracy**: Percentage of correct predictions
- **Precision**: True positives / (True positives + False positives)
- **Recall**: True positives / (True positives + False negatives)
- **F1-score**: Harmonic mean of precision and recall
- **Cross-validation**: 3-fold stratified CV for robust estimation

## Performance Tips

### 1. **Optimal k Selection**
- **Small k (1-3)**: More sensitive to noise, higher variance
- **Medium k (5-9)**: Often optimal balance
- **Large k (11+)**: Smoother decision boundaries, may underfit

### 2. **Data Quality**
- Ensure distance matrices are symmetric
- Check for missing annotations
- Verify protein ID consistency between matrices and annotations

### 3. **Memory Considerations**
- Large distance matrices (>1000×1000) may require significant RAM
- Consider batch processing for very large datasets

## Troubleshooting

### Common Issues

#### 1. **Import Errors**
```bash
# Install missing packages
pip install scikit-learn pandas numpy matplotlib seaborn biopython
```

#### 2. **File Not Found**
```bash
# Check file paths
ls -la ./fastas/
ls -la *.csv
```

#### 3. **Memory Issues**
```python
# For large matrices, consider:
# - Reducing the number of proteins
# - Using float32 instead of float64
# - Processing in smaller batches
```

#### 4. **Poor Performance**
- Check class balance in your dataset
- Verify distance matrix quality
- Try different k values
- Consider feature selection or dimensionality reduction

### Getting Help

1. **Check the output logs** for detailed error messages
2. **Verify input file formats** match the expected structure
3. **Test with smaller datasets** first
4. **Check EC level availability** in your annotation files

## Citation

If you use this k-NN classifier in your research, please cite:

```bibtex
@software{protsec_knn_classifier,
  title={k-NN Classifier for Enzyme Class Prediction using ProtSEC Distance Matrices},
  author={Your Name},
  year={2025},
  url={https://github.com/omics-lab/ProtSEC}
}
```

## License

This software is available under the same license as the ProtSEC project.
