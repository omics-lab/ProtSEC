library(umap)
library(tidyverse)

# Check current working directory
setwd("//wsl$/Ubuntu/home/rashedul/project/ProSEC/data/clustering")
getwd()

# ---- UMAP phosphatase----

# UMAP for multiple methods

# List all CSV files in the protein_kinase directory
csv_files <- list.files("protein_kinase", pattern = "\\_matrix\\.csv$", full.names = TRUE)

# Extract method names from file names
method_names <- gsub("protein_kinase/|kinase\\.|_matrix\\.csv", "", basename(csv_files))

# Initialize empty list to store UMAP dataframes
umap_list <- list()

# Process each CSV file
for (i in seq_along(csv_files)) {
  # Read distance matrix
  dist_mat <- as.matrix(read.csv(csv_files[i], row.names = 1))
  
  # Perform UMAP
  umap_result <- umap(dist_mat)
  
  # Create dataframe
  umap_df <- as.data.frame(umap_result$layout)
  umap_df$accession <- rownames(dist_mat)
  umap_df$method <- method_names[i]  # Add method column
  
  # Store in list
  umap_list[[i]] <- umap_df
}

# Combine all UMAP dataframes
combined_umap_df <- do.call(rbind, umap_list)

# Add labels and filtering
labels_df <- read.csv("protein_kinase/kinase_labels.csv", stringsAsFactors = FALSE)
combined_umap_labeled <- merge(combined_umap_df, labels_df, by = "accession")

# Apply exclusions and filtering
exclude_df <- read.table("protein_kinase/exclude.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(exclude_df) <- c("accession", "label")
combined_umap_labeled <- subset(combined_umap_labeled, !(accession %in% exclude_df$accession)) %>%
  filter(label != "")

# Create combined plot
p_combined <- ggplot(combined_umap_labeled, aes(x = V1, y = V2, color = label)) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Paired") +
  labs(title = "UMAP by Method", x = "UMAP 1", y = "UMAP 2", color = "Label") +
  facet_wrap(~method, ncol = 3) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "gray", size = 1),
        strip.background = element_rect(fill="white"),
        legend.position = "bottom",
        axis.text = element_text(color = "black", vjust = 0.5, hjust=1))

ggsave("../plots/combined_kinase_umap.pdf", plot = p_combined, width = 15, height = 10)

# Calculate trustworthiness for all methods
trustworthiness_results <- data.frame(
  method = character(),
  trustworthiness = numeric(),
  stringsAsFactors = FALSE
)

# Process each CSV file and calculate trustworthiness
for (i in seq_along(csv_files)) {
  # Read distance matrix
  dist_mat <- as.matrix(read.csv(csv_files[i], row.names = 1))
  
  # Perform UMAP
  umap_result <- umap(dist_mat)
  
  # Calculate trustworthiness
  X <- as.matrix(dist_mat) 
  Y <- umap_result$layout 
  trust_score <- calculate_trustworthiness(X, Y, k = 15)
  
  # Store results
  trustworthiness_results <- rbind(trustworthiness_results, 
                                   data.frame(method = method_names[i], 
                                             trustworthiness = trust_score))
  
  cat("Method:", method_names[i], "- Trustworthiness:", round(trust_score, 3), "\n")
}

# Print summary table
print("Trustworthiness Summary:")
print(trustworthiness_results)

# Save results
write.csv(trustworthiness_results, "../protein_kinase/trustworthiness_kinase_results.csv", row.names = FALSE)

# Create trustworthiness comparison plot
p_trust <- ggplot(trustworthiness_results, aes(x = reorder(method, trustworthiness), y = trustworthiness)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = round(trustworthiness, 3)), vjust = -0.5) +
  theme_minimal() +
  labs(title = "Trustworthiness Comparison Across Methods",
       x = "Method", 
       y = "Trustworthiness Score",
       subtitle = "Higher scores indicate better preservation of local neighborhoods") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

ggsave("../plots/trustworthiness_kinase_comparison.pdf", plot = p_trust, width = 10, height = 6)


# ---- UMAP for single method (ProtSEC) ----

# single file 
dist_mat <- as.matrix(read.csv("protein_kinase/kinase.proseq_dist_matrix.csv", row.names = 1))
umap_result <- umap(dist_mat)
umap_df <- as.data.frame(umap_result$layout)
umap_df$accession <- rownames(corr_mat)

labels_df <- read.csv("protein_kinase/kinase_labels.csv", stringsAsFactors = FALSE)
umap_df_labeled <- merge(umap_df, labels_df, by = "accession")

exclude_df <- read.table("protein_kinase/exclude.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(exclude_df) <- c("accession", "label")
umap_df_labeled <- subset(umap_df_labeled, !(accession %in% exclude_df$accession)) %>% filter(label != "")

p_kinase <- ggplot(umap_df_labeled, aes(x = V1, y = V2, color = label)) +
  geom_point(size = 1.5) +
  # geom_text(vjust = 1.5, hjust = 1.1, size = 2.5) +
  scale_color_brewer(palette = "Paired") +
  labs(title = "ProtSEC", x = "UMAP 1", y = "UMAP 2", color = "Label") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "gray", size = 1),
        strip.background =element_rect(fill="white"),
        legend.position = "bottom",
        axis.text = element_text(color = "black", vjust = 0.5, hjust=1))
ggsave("../plots/kinase_umap.pdf", plot = p_kinase, width = 7, height = 7)


# trustworthiness between original distance and UMAP layout
X <- as.matrix(dist_mat) 
Y <- umap_result$layout 
calculate_trustworthiness(X, Y, k = 15)


calculate_trustworthiness(X, Y, k = 15)
# 0.90 for phosphatase
# 0.91 for kinase
# 0.67 for sam

# trustworthiness function
calculate_trustworthiness <- function(X, Y, k = 15) {
  library(FNN)

  n <- nrow(X)

  # Get k-NN from original and embedded space
  nn_X <- get.knn(X, k = n - 1)$nn.index # full ranking
  nn_Y <- get.knn(Y, k = k)$nn.index # just top k

  t_sum <- 0

  for (i in 1:n) {
    # Actual ranks of neighbors in original space
    true_neighbors <- nn_X[i, ]
    y_neighbors <- nn_Y[i, ]

    for (j in y_neighbors) {
      rank_j <- which(true_neighbors == j)
      if (length(rank_j) == 0) {
        next # skip if not found
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


