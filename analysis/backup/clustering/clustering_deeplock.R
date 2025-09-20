library(umap)
library(tidyverse)
library(patchwork)

# Check current working directory
setwd("//wsl$/Ubuntu/home/rashedul/project/ProSEC/data/clustering")
getwd()

# cd /home/rashedul/project/ProSEC
# python3 get_phase_dist_mat.py -i data/clustering/deeploc/deeploc_data.fasta -o  data/clustering/deeploc/deeploc_data.proseq_matrix.csv

# ---- UMAP deeploc----

# List all CSV files in the Deeploc directory
csv_files <- list.files("deeploc", pattern = "\\_matrix\\.csv$", full.names = TRUE)

# Extract method names from file names
method_names <- gsub("deeploc/|deeploc\\.|_matrix\\.csv", "", basename(csv_files))

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
dim(combined_umap_df)
head(combined_umap_df)
head(combined_umap_df$accession)

# Read sequence descriptions from FASTA file
fasta_lines <- readLines("deeploc/deeploc_data.fasta")
header_lines <- fasta_lines[grepl("^>", fasta_lines)]

# Extract accession and location from FASTA headers
# Remove "test" from header lines first
clean_header_lines <- gsub(" test", "", header_lines)

# Split by space to get accession and location
header_parts <- strsplit(clean_header_lines, " ")
labels_df <- data.frame(
  accession = sapply(header_parts, function(x) sub("^>", "", x[1])),
  location = sapply(header_parts, function(x) if(length(x) > 1) x[2] else ""),
  stringsAsFactors = FALSE
)

# Remove characters after dash (e.g., '-S', '-M', '-U')
labels_df$location <- sub("-.*$", "", labels_df$location)
unique(labels_df$location)
head(labels_df)

combined_umap_labeled <- merge(combined_umap_df, labels_df, by = "accession")
head(combined_umap_labeled)

# Filter out values outside the range -12 to 12 for V1 and V2
combined_umap_labeled <- combined_umap_labeled %>%
  filter(V1 >= -12 & V1 <= 12 & V2 >= -12 & V2 <= 12) # protseq generates some large values

# Order methods by specified order
method_order <- c("ProtSEC", "esm2_large", "esm2_small", "prot_t5", "prot_bert")
combined_umap_labeled$method <- factor(combined_umap_labeled$method, levels = method_order)

# Create combined plot
p_combined <- ggplot(combined_umap_labeled, aes(x = V1, y = V2, color = location)) +
  geom_point(size = 0.3, alpha = 0.6) +  # Smaller points with transparency
  scale_color_brewer(palette = "Paired") +  
  labs(title = "DeepLoc", x = "UMAP 1", y = "UMAP 2", color = "Location") +
  facet_wrap(~method, ncol = 3, scales = "free") +
  theme_minimal() +  
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "gray", size = 0.5),
        strip.background = element_rect(fill="white"),
        legend.position = "bottom",
        axis.text = element_text(color = "black", size = 8),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))  

ggsave("../plots/combined_deeploc_umap.pdf", plot = p_combined, width = 15, height = 10, dpi = 150)

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

# Create trustworthiness comparison plot
trustworthiness_results$method <- factor(trustworthiness_results$method, levels = method_order)

p_trust <- ggplot(trustworthiness_results, aes(x = reorder(method, trustworthiness), y = trustworthiness)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = round(trustworthiness, 3)), vjust = -0.5) +
  labs(x = "", 
       y = "Trustworthiness Score") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "gray", size = 1),
        strip.background = element_rect(fill="white"),
        legend.position = "bottom",
        axis.text = element_text(color = "black", angle = 45, vjust = 0.5, hjust=1)) +
  ylim(0, 1.1)

ggsave("../plots/trustworthiness_deeploc_comparison.pdf", plot = p_trust, width = 10, height = 6)

# Create the combined layout using patchwork
combined_layout <- p_combined + 
  inset_element(p_trust, 
                left = 0.68, bottom = -0.11, right = 1, top = 0.49)

ggsave("../plots/combined_deeploc_umap_with_trust.pdf", plot = combined_layout, width = 15, height = 10)


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
