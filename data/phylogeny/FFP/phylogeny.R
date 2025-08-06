# .libPaths("C:/Program Files/R/R-4.3.1/library")

library(tidyverse)
library(ape)  
library(phangorn)   
#BiocManager::install("msa") 
library(msa)
library(Biostrings) 

# Check current working directory
setwd("//wsl$/Ubuntu/home/rashedul/project/ProSEC/data/phylogeny")


# ---- ClustalW + UPGMA for two sequence sets ----
seqs1 <- readDNAStringSet("./FFP/10-BetaSet_edited.fasta")
alignment1 <- msa(seqs1, method = "ClustalW")
alignment_phangorn1 <- msaConvert(alignment1, type = "phangorn::phyDat")
dist_matrix1 <- dist.ml(alignment_phangorn1)
tree1 <- upgma(dist_matrix1)

seqs2 <- readDNAStringSet("./FFP/17-BetaSet_edited.fasta")
alignment2 <- msa(seqs2, method = "ClustalW")
alignment_phangorn2 <- msaConvert(alignment2, type = "phangorn::phyDat")
dist_matrix2 <- dist.ml(alignment_phangorn2)
tree2 <- upgma(dist_matrix2)

pdf("../plots/clustalw_upgma_twoSeqs.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
plot(tree1, main = "10-BetaSet + UPGMA", cex = 0.8)
plot(tree2, main = "17-BetaSet + UPGMA", cex = 0.8)
dev.off()


# ---- Bootstrap Analysis ----
# did not work, need to fix later


# ---- ProSEQ + UPGMA ----
generate_upgma_tree <- function(matrix_csv, method_name) {
  corr_mat <- as.matrix(read.csv(matrix_csv, row.names = 1))
  dist_mat <- as.dist(1 - corr_mat)
  tree <- upgma(dist_mat)
  assign(paste0("tree_", method_name, "_upgma"), tree, envir = .GlobalEnv)
  return(tree)
}

# ---- Plotting the trees ----
matrix_files <- list.files("./FFP", pattern = "proseq_matrix.csv$", full.names = TRUE)
method_names <- gsub("10-BetaSet.|\\.csv", "", basename(matrix_files))

pdf("../plots/Protseq_2BetaSet_UPGMA_Trees.pdf", width = 16, height = 8)
par(mfrow = c(1, 2)) 

for (i in seq_along(matrix_files)) {
  tree <- generate_upgma_tree(matrix_files[i], method_names[i])
  plot(tree, main = paste(method_names[i], "+ UPGMA"), cex = 0.8)
}
dev.off()


# ---- PLMs + UPGMA ----
# csvs are distance matrices already
generate_upgma_tree <- function(matrix_csv, method_name) {
  corr_mat <- as.matrix(read.csv(matrix_csv, row.names = 1))
  # dist_mat <- as.dist(1 - corr_mat)
  tree <- upgma(corr_mat)
  assign(paste0("tree_", method_name, "_upgma"), tree, envir = .GlobalEnv)
  return(tree)
}

# ---- Plotting the trees ----
matrix_files <- list.files("./FFP", pattern = "_matrix.csv$", full.names = TRUE)
matrix_files <- matrix_files[!grepl("proseq_matrix.csv", matrix_files)]
method_names <- gsub("10-BetaSet.|\\.csv", "", basename(matrix_files))

pdf("../plots/all_4PLMs_BetaSet_UPGMA_Trees.pdf", width = 16, height = 8)
par(mfrow = c(2, 2)) 

for (i in seq_along(matrix_files)) {
  tree <- generate_upgma_tree(matrix_files[i], method_names[i])
  plot(tree, main = paste(method_names[i], "+ UPGMA"), cex = 0.8)
}
dev.off()


# tree_proseq_upgma <- generate_upgma_tree("./FFP/10-BetaSet.proseq_matrix.csv", "proseq")
test <- generate_upgma_tree("./FFP/esm2_small_cosine_dis_matrix.csv", "esm2_small")
pdf("../plots/test.pdf", width = 8, height = 6)
plot(test, main = "esm2_small + UPGMA", cex = 0.8)
dev.off()




