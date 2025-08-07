library(tidyverse)
library(ape)  
library(phangorn)   
#BiocManager::install("msa") 
library(msa)
library(Biostrings)
BiocManager::install("ggtree") 
library(ggtree)
library(patchwork) 


# Check current working directory
setwd("//wsl$/Ubuntu/home/rashedul/project/ProSEC/data/phylogeny")


# ---- ClustalW + UPGMA for two sequence sets ----
seqs1 <- readDNAStringSet("./FFP/10-BetaSet_edited.fasta")
alignment1 <- msa(seqs1, method = "ClustalW")
alignment_phangorn1 <- msaConvert(alignment1, type = "phangorn::phyDat")
dist_matrix1 <- dist.ml(alignment_phangorn1)
tree_clustalw_upgma_1 <- upgma(dist_matrix1)

seqs2 <- readDNAStringSet("./FFP/17-BetaSet_edited.fasta")
alignment2 <- msa(seqs2, method = "ClustalW")
alignment_phangorn2 <- msaConvert(alignment2, type = "phangorn::phyDat")
dist_matrix2 <- dist.ml(alignment_phangorn2)
tree_clustalw_upgma_2 <- upgma(dist_matrix2)

# Create ggtree plot objects for clustalW trees
p1 <- ggtree(tree_clustalw_upgma_1) + 
  geom_tiplab(size = 4, hjust = -0.1) + 
  ggtitle("ClustalW") +
  theme_tree2() +
  xlim(0, max(dist.nodes(tree_clustalw_upgma_1)) * 0.65)

p2 <- ggtree(tree_clustalw_upgma_2) + 
  geom_tiplab(size = 4, hjust = -0.1) + 
  ggtitle("ClustalW") +
  theme_tree2() +
  xlim(0, max(dist.nodes(tree_clustalw_upgma_2)) * 0.65)

# ---- PLMs + ProtSEC plotting ----

generate_upgma_tree_dist <- function(matrix_csv, method_name) {
  dist_mat <- as.matrix(read.csv(matrix_csv, row.names = 1))
  tree <- upgma(dist_mat)
  assign(paste0("tree_", method_name, "_upgma"), tree, envir = .GlobalEnv)
  return(tree)
}

# Plotting the trees for 10-BetaSet
matrix_files <- list.files("./FFP", pattern = "10-BetaSet.*_matrix.csv$", full.names = TRUE)
method_names <- gsub("10-BetaSet\\.|17-BetaSet\\.|_matrix\\.csv", "", basename(matrix_files))

# Create a list to store ggtree plots
plot_list <- list()

# Define the desired order for the plots
desired_order <- c("ClustalW", "protseq", "esm2_large", "esm2_small", "prot_t5", "prot_bert")

# Add ClustalW (p1) as the first plot
plot_list[["ClustalW"]] <- p1

# Generate plots for each method in the desired order
for (method in desired_order[-1]) {  # Skip "ClustalW" since it's already added
  # Find the corresponding matrix file for this method
  matrix_file <- matrix_files[grepl(paste0(method, "_matrix"), matrix_files)]
  
  if (length(matrix_file) > 0) {
    tree <- generate_upgma_tree_dist(matrix_file[1], method)
    plot_list[[method]] <- ggtree(tree) + 
      geom_tiplab(size = 4, hjust = -0.1) + 
      ggtitle(paste(method)) +
      theme_tree2() +
      theme(plot.title = element_text()) +
      xlim(0, max(dist.nodes(tree)) * 0.65)
  }
}

# Convert named list to ordered list based on desired_order
plot_list_ordered <- plot_list[desired_order[desired_order %in% names(plot_list)]]

# Arrange plots in a 2x3 grid and save directly to PDF
combined_plot <- wrap_plots(plot_list_ordered, ncol = 3, nrow = 2)
ggsave("../plots/all_4PLMs_10BetaSet_UPGMA_Trees.pdf", plot = combined_plot, width = 20, height = 10)


## Plotting the trees for 17-BetaSet
matrix_files <- list.files("./FFP", pattern = "17-BetaSet.*_matrix.csv$", full.names = TRUE)
method_names <- gsub("10-BetaSet\\.|17-BetaSet\\.|_matrix\\.csv", "", basename(matrix_files))

# Create a list to store ggtree plots
plot_list <- list()

# Define the desired order for the plots
desired_order <- c("ClustalW", "protseq", "esm2_large", "esm2_small", "prot_t5", "prot_bert")

# Add ClustalW (p1) as the first plot
plot_list[["ClustalW"]] <- p2

# Generate plots for each method in the desired order
for (method in desired_order[-1]) {  # Skip "ClustalW" since it's already added
  # Find the corresponding matrix file for this method
  matrix_file <- matrix_files[grepl(paste0(method, "_matrix"), matrix_files)]
  
  if (length(matrix_file) > 0) {
    tree <- generate_upgma_tree_dist(matrix_file[1], method)
    plot_list[[method]] <- ggtree(tree) + 
      geom_tiplab(size = 4, hjust = -0.1) + 
      ggtitle(paste(method)) +
      theme_tree2() +
      theme(plot.title = element_text()) +
      xlim(0, max(dist.nodes(tree)) * 0.65)
  }
}

# Convert named list to ordered list based on desired_order
plot_list_ordered <- plot_list[desired_order[desired_order %in% names(plot_list)]]

# Arrange plots in a 2x3 grid and save directly to PDF
combined_plot <- wrap_plots(plot_list_ordered, ncol = 3, nrow = 2)
ggsave("../plots/all_4PLMs_17BetaSet_UPGMA_Trees.pdf", plot = combined_plot, width = 20, height = 10)


# ---- Branch Score Distance (BSD) analysis ----

# Function to normalize tree edge lengths
normalize_tree <- function(tree) {
  tree$edge.length <- tree$edge.length / sum(tree$edge.length)
  return(tree)
}

# List matrix files for each BetaSet 
matrix_files_10 <- list.files("./FFP", pattern = "10-BetaSet.*_matrix.csv$", full.names = TRUE)

matrix_files_17 <- list.files("./FFP", pattern = "17-BetaSet.*_matrix.csv$", full.names = TRUE)

method_names_10 <- gsub("10-BetaSet.|\\.csv", "", basename(matrix_files_10))
method_names_17 <- gsub("17-BetaSet.|\\.csv", "", basename(matrix_files_17))

# 10-BetaSet BSD collection
bsd_10 <- numeric(length(matrix_files_10))
for (i in seq_along(matrix_files_10)) {
  tree <- generate_upgma_tree_dist(matrix_files_10[i], method_names_10[i])
  tree_norm <- normalize_tree(tree)
  ref_norm <- normalize_tree(tree_clustalw_upgma_1)
  bsd_10[i] <- dist.topo(ref_norm, tree_norm, method = "score")
  cat(method_names_10[i], "Normalized BSD:", bsd_10[i], "\n")
}

bsd_10_df <- data.frame(Method = method_names_10, Normalized_BSD = bsd_10, Source = "bsd_10")

# 17-BetaSet BSD collection
bsd_17 <- numeric(length(matrix_files_17))
for (i in seq_along(matrix_files_17)) {
  tree <- generate_upgma_tree_dist(matrix_files_17[i], method_names_17[i])
  tree_norm <- normalize_tree(tree)
  ref_norm <- normalize_tree(tree_clustalw_upgma_2)
  bsd_17[i] <- dist.topo(ref_norm, tree_norm, method = "score")
  cat(method_names_17[i], "Normalized BSD:", bsd_17[i], "\n")
}

bsd_17_df <- data.frame(Method = method_names_17, Normalized_BSD = bsd_17, Source = "bsd_17")

all_bsd <- rbind(bsd_10_df, bsd_17_df)
all_bsd$Method <- gsub("_matrix", "", all_bsd$Method)
write.csv(all_bsd, "./FFP/all_bsd.csv", row.names = FALSE)

# Define the desired order for Methods (updated without "_matrix")
method_order <- c("protseq", "esm2_large", "esm2_small", "prot_t5", "prot_bert")

# Convert Method to factor with specified order
all_bsd$Method <- factor(all_bsd$Method, levels = method_order)

# Barplot for bsd_10
p_10 <- ggplot(subset(all_bsd, Source == "bsd_10"), aes(x = Method, y = Normalized_BSD, fill = Method)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  labs(
    x = "Method",
    y = "Normalized BSD",
    fill = "Method",
    title = "Normalized BSD for 10-BetaSet"
  )

ggsave("../plots/bsd_10_barplot.pdf", plot = p_10, width = 8, height = 8)

# Barplot for bsd_17
p_17 <- ggplot(subset(all_bsd, Source == "bsd_17"), aes(x = Method, y = Normalized_BSD, fill = Method)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  labs(
    x = "Method",
    y = "Normalized_BSD",
    fill = "Method",
    title = "Normalized BSD for 17-BetaSet"
  )

ggsave("../plots/bsd_17_barplot.pdf", plot = p_17, width = 8, height = 8)

# Calculate average BSD by Source excluding protseq
avg_bsd_by_source <- all_bsd %>%
  filter(Method != "protseq") %>%
  group_by(Source) %>%
  summarise(Average_Normalized_BSD = mean(Normalized_BSD),
            .groups = "drop")

print("Average Normalized BSD by Source (excluding protseq):")
print(avg_bsd_by_source)


# ---- Backup Analysis ----


# Bootstrap 
# did not work, need to fix later


# # ---- ProSEQ corr converted to dist ----
# generate_upgma_tree_corr <- function(matrix_csv, method_name) {
#   corr_mat <- as.matrix(read.csv(matrix_csv, row.names = 1))
#   dist_mat <- as.dist(1 - corr_mat)
  
#   # Save distance matrix as CSV
#   write.csv(as.matrix(dist_mat), paste0("dist_", method_name, ".csv"))
  
#   tree <- upgma(dist_mat)
#   assign(paste0("tree_", method_name, "_upgma"), tree, envir = .GlobalEnv)
#   return(tree)
# }

# # ---- Plotting the trees ----
# matrix_files <- list.files("./FFP/temp_corr/", pattern = "proseq_matrix.csv$", full.names = TRUE)
# method_names <- gsub("10-BetaSet.|\\.csv", "", basename(matrix_files))

# pdf("../plots/Protseq_2BetaSet_UPGMA_Trees.pdf", width = 16, height = 8)
# par(mfrow = c(1, 2)) 

# for (i in seq_along(matrix_files)) {
#   tree <- generate_upgma_tree_corr(matrix_files[i], method_names[i])
#   plot(tree, main = paste(method_names[i], "+ UPGMA"), cex = 0.8)
# }
# dev.off()


# # List ProtSeq matrix files for each BetaSet
# matrix_files_10 <- list.files("./FFP", pattern = "10-BetaSet.*proseq_matrix.csv$", full.names = TRUE)
# matrix_files_17 <- list.files("./FFP", pattern = "17-BetaSet.*proseq_matrix.csv$", full.names = TRUE)

# method_names_10 <- gsub("10-BetaSet.|\\.csv", "", basename(matrix_files_10))
# method_names_17 <- gsub("17-BetaSet.|\\.csv", "", basename(matrix_files_17))

# # 10-BetaSet with normalization-
# cat("Normalized Branch Score Distances for 10-BetaSet:\n")
# for (i in seq_along(matrix_files_10)) {
#   tree <- generate_upgma_tree_corr(matrix_files_10[i], method_names_10[i])
#   tree_norm <- normalize_tree(tree)
#   ref_norm <- normalize_tree(tree_clustalw_upgma_1)
#   bsd <- dist.topo(ref_norm, tree_norm, method = "score")
#   cat(method_names_10[i], "Normalized BSD:", bsd, "\n")
# }

# protseq_bsd_10_df <- data.frame(Method = method_names_10, Normalized_BSD = bsd, Source = "bsd_10")

# # 17-BetaSet with normalization
# cat("\nNormalized Branch Score Distances for 17-BetaSet:\n")
# for (i in seq_along(matrix_files_17)) {
#   tree <- generate_upgma_tree_corr(matrix_files_17[i], method_names_17[i])
#   tree_norm <- normalize_tree(tree)
#   ref_norm <- normalize_tree(tree_clustalw_upgma_2)
#   bsd <- dist.topo(ref_norm, tree_norm, method = "score")
#   cat(method_names_17[i], "Normalized BSD:", bsd, "\n")
# }

# protseq_bsd_17_df <- data.frame(Method = method_names_10, Normalized_BSD = bsd, Source = "bsd_17")
