## knn

setwd("/Users/rashedulislam/Documents/git_repos/ProtSEC/analysis/EC_CARE_task1/")

dist_mat = read.csv("./fastas/sample_2000_ProtSEC_matrix.csv",  row.names = 1)
print(dist_mat[1:5, 1:5])
dim(dist_mat)

## Prepare training data
training_metadata = read.csv("protein_train.csv")
print(training_metadata[1:5, 1:5])

# Get intersection of distance matrix and training data entries
dist_mat_proteins <- rownames(dist_mat)
training_proteins <- training_metadata$Entry

# Find common proteins
common_proteins <- intersect(dist_mat_proteins, training_proteins)
length(common_proteins)

# Filter distance matrix to keep only rows for common proteins (all columns)
train_mat <- dist_mat[common_proteins, ]
print(train_mat[1:5, 1:5]) 
dim(train_mat)

# Extract labels in the same order as train_mat rows
train_labels <- training_metadata$EC1[match(rownames(train_mat), training_metadata$Entry)]
length(train_labels)

## Prepare test data 
test_metadata = read.csv("30_protein_test.csv")
print(test_metadata[1:5, 1:5])
dim(test_metadata) 

# Get intersection of distance matrix and test data entries
test_proteins <- test_metadata$Entry

# Find common proteins between test data and distance matrix
common_test_proteins <- intersect(dist_mat_proteins, test_proteins)
length(common_test_proteins)

# Filter distance matrix to keep only rows for common test proteins (all columns)
test_mat <- dist_mat[common_test_proteins, ]
print(test_mat[1:5, 1:5]) 
dim(test_mat)

# Extract test labels in the same order as test_mat rows
test_labels <- test_metadata$EC1[match(rownames(test_mat), test_metadata$Entry)]
length(test_labels)
print("Test labels length:")
print(length(test_labels))
print("Test matrix rows:")
print(nrow(test_mat))

# Run k-NN
knn_pred <- knn(train = train_mat, test = test_mat, cl = train_labels, k = 5)
print("k-NN predictions:")
print(knn_pred)
print("Actual labels:")
print(test_labels)

# Prediction testing codes
library(caret)

# Create confusion matrix
conf_matrix <- table(Predicted = knn_pred, Actual = test_labels)
print(conf_matrix)

# Calculate accuracy
total_predictions <- sum(conf_matrix)
correct_predictions <- sum(diag(conf_matrix))
accuracy <- (correct_predictions / total_predictions) * 100

print(paste("Total predictions:", total_predictions))
print(paste("Correct predictions:", correct_predictions))
print(paste("Accuracy:", round(accuracy, 2), "%"))


## Example code
# Create example matrix for k-NN testing
library(class)

# Create a simple 5x5 example distance matrix
example_dist <- matrix(c(
  0.0, 0.2, 0.8, 0.9, 0.7,
  0.2, 0.0, 0.6, 0.8, 0.5,
  0.8, 0.6, 0.0, 0.3, 0.4,
  0.9, 0.8, 0.3, 0.0, 0.2,
  0.7, 0.5, 0.4, 0.2, 0.0
), nrow = 5, byrow = TRUE)

# Add row and column names
rownames(example_dist) <- c("P1", "P2", "P3", "P4", "P5")
colnames(example_dist) <- c("P1", "P2", "P3", "P4", "P5")
print(example_dist) 

# Example labels
example_labels <- c("A", "A", "B", "B", "B")

# Split into train (first 3) and test (last 2)
train_indices <- 1:3
test_indices <- 4:5

X_train <- example_dist[train_indices, train_indices]
X_test <- example_dist[test_indices, train_indices]
y_train <- example_labels[train_indices]
y_test <- example_labels[test_indices]

print(X_train)
print(X_test)
print(y_train)
print(y_test)

# Run k-NN
knn_pred <- knn(train = X_train, test = X_test, cl = y_train, k = 1)
print("k-NN predictions:")
print(knn_pred)
print("Actual labels:")
print(y_test)

# Prediction testing codes
library(caret)

# Create confusion matrix
conf_matrix <- table(Predicted = knn_pred, Actual = y_test)
conf_matrix
