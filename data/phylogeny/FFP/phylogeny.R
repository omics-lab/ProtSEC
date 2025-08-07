# .libPaths("C:/Program Files/R/R-4.3.1/library")

library(tidyverse)
library(ape)  
library(phangorn)   
#BiocManager::install("msa") 
library(msa)
library(Biostrings) 
library(ggplot2)

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

pdf("../plots/clustalw_upgma_twoSeqs.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
plot(tree_clustalw_upgma_1, main = "10-BetaSet + UPGMA", cex = 0.8)
plot(tree_clustalw_upgma_2, main = "17-BetaSet + UPGMA", cex = 0.8)
dev.off()


# ---- PLMs + ProtSEC plotting ----
# csvs are distance matrices already
generate_upgma_tree_dist <- function(matrix_csv, method_name) {
  dist_mat <- as.matrix(read.csv(matrix_csv, row.names = 1))
  tree <- upgma(dist_mat)
  assign(paste0("tree_", method_name, "_upgma"), tree, envir = .GlobalEnv)
  return(tree)
}

# Plotting the trees 
matrix_files <- list.files("./FFP", pattern = "_matrix.csv$", full.names = TRUE)
method_names <- gsub("10-BetaSet\\.|17-BetaSet\\.|_matrix\\.csv", "", basename(matrix_files))

pdf("../plots/all_4PLMs_BetaSet_UPGMA_Trees.pdf", width = 16, height = 8)
par(mfrow = c(2, 3)) 

for (i in seq_along(matrix_files)) {
  tree <- generate_upgma_tree_dist(matrix_files[i], method_names[i])
  plot(tree, main = paste(method_names[i], "+ UPGMA"), cex = 0.8)
}
dev.off()


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
  labs(
    x = "Method",
    y = "Normalized BSD",
    fill = "Method",
    title = "Normalized BSD for 10-BetaSet"
  )

ggsave("../plots/bsd_10_barplot.pdf", plot = p_10, width = 8, height = 5)

# Barplot for bsd_17
p_17 <- ggplot(subset(all_bsd, Source == "bsd_17"), aes(x = Method, y = Normalized_BSD, fill = Method)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Method",
    y = "Normalized BSD",
    fill = "Method",
    title = "Normalized BSD for 17-BetaSet"
  )

ggsave("../plots/bsd_17_barplot.pdf", plot = p_17, width = 8, height = 5)


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
