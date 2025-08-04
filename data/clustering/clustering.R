library(umap)
library(tidyverse)

# Check current working directory
setwd("//wsl$/Ubuntu/home/rashedul/project/ProSEC/data/clustering")

# ---- UMAP phosphatase----

corr_mat <- as.matrix(read.csv("phosphatase/phosphatase.proseq_matrix.csv", row.names = 1))
dist_mat <- 1 - corr_mat

umap_result <- umap(dist_mat)
umap_df <- as.data.frame(umap_result$layout)
umap_df$accession <- rownames(corr_mat)

labels_df <- read.csv("phosphatase/phosphatase_labels.csv", stringsAsFactors = FALSE)
labels_df$accession <- gsub("\\|.*", "", labels_df$accession)
umap_df$accession <- gsub("\\|.*", "", umap_df$accession)
umap_df_labeled <- merge(umap_df, labels_df, by = "accession")

exclude_df <- read.table("phosphatase/exclude.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(exclude_df) <- c("accession", "label")
exclude_df$accession <- gsub("\\|.*", "", exclude_df$accession)
umap_df_labeled <- subset(umap_df_labeled, !(accession %in% exclude_df$accession))

ggplot(umap_df_labeled, aes(x = V1, y = V2, color = label, label = label)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 1.1, size = 2.5) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  labs(title = "UMAP Colored by Label", x = "UMAP 1", y = "UMAP 2", color = "Label")

# Example: trustworthiness between original distance and UMAP layout
X <- as.matrix(dist_mat)         # original distance or features
Y <- umap_result$layout          # UMAP layout

calculate_trustworthiness(X, Y, k = 15)

# ---- UMAP kinase----

corr_mat <- as.matrix(read.csv("protein_kinase/kinase.proseq_matrix.csv", row.names = 1))
dist_mat <- 1 - corr_mat

umap_result <- umap(dist_mat)
umap_df <- as.data.frame(umap_result$layout)
umap_df$accession <- rownames(corr_mat)

labels_df <- read.csv("protein_kinase/kinase_labels.csv", stringsAsFactors = FALSE)
labels_df$accession <- gsub("\\|.*", "", labels_df$accession)
umap_df$accession <- gsub("\\|.*", "", umap_df$accession)
umap_df_labeled <- merge(umap_df, labels_df, by = "accession")

exclude_df <- read.table("protein_kinase/exclude.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(exclude_df) <- c("accession", "label")
exclude_df$accession <- gsub("\\|.*", "", exclude_df$accession)
umap_df_labeled <- subset(umap_df_labeled, !(accession %in% exclude_df$accession))

ggplot(umap_df_labeled, aes(x = V1, y = V2, color = label, label = label)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 1.1, size = 2.5) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  labs(title = "UMAP Colored by Label", x = "UMAP 1", y = "UMAP 2", color = "Label")

# Example: trustworthiness between original distance and UMAP layout
X <- as.matrix(dist_mat)         # original distance or features
Y <- umap_result$layout          # UMAP layout

calculate_trustworthiness(X, Y, k = 15)

# ---- UMAP sam----

corr_mat <- as.matrix(read.csv("radical_sam/radicalsam.proseq_matrix.csv", row.names = 1))
dist_mat <- 1 - corr_mat
umap_result <- umap(dist_mat)
umap_df <- as.data.frame(umap_result$layout)
umap_df$accession <- rownames(corr_mat)

labels_df <- read.csv("radical_sam/radicalsam_labels.csv", stringsAsFactors = FALSE)
labels_df$accession <- gsub("\\|.*", "", labels_df$accession)
umap_df$accession <- gsub("\\|.*", "", umap_df$accession)

umap_df_labeled <- merge(umap_df, labels_df, by = "accession")
head(umap_df_labeled)

ggplot(umap_df_labeled, aes(x = V1, y = V2, color = label, label = label)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 1.1, size = 2.5) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  labs(title = "UMAP of Phosphatases Colored by Label", x = "UMAP 1", y = "UMAP 2", color = "Label")


# Example: trustworthiness between original distance and UMAP layout
X <- as.matrix(dist_mat)         # original distance or features
Y <- umap_result$layout          # UMAP layout

calculate_trustworthiness(X, Y, k = 15)
# 0.90 for phosphatase
# 0.91 for kinase
# 0.67 for sam

# trustworthiness function 
calculate_trustworthiness <- function(X, Y, k = 15) {
  library(FNN)
  
  n <- nrow(X)
  
  # Get k-NN from original and embedded space
  nn_X <- get.knn(X, k = n - 1)$nn.index  # full ranking
  nn_Y <- get.knn(Y, k = k)$nn.index      # just top k
  
  t_sum <- 0
  
  for (i in 1:n) {
    # Actual ranks of neighbors in original space
    true_neighbors <- nn_X[i, ]
    y_neighbors <- nn_Y[i, ]
    
    for (j in y_neighbors) {
      rank_j <- which(true_neighbors == j)
      if (length(rank_j) == 0) {
        next  # skip if not found
      }
      if (rank_j > k) {
        t_sum <- t_sum + (rank_j - k)
      }
    }
  }
  
  # Trustworthiness formula
  t <- 1 - (2 / (n * k * (2 * n - 3 * k - 1))) * t_sum
  return(t)
}
