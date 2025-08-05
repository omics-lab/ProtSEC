# .libPaths("C:/Program Files/R/R-4.3.1/library")

library(tidyverse)
library(ape)  
library(phangorn)   
#BiocManager::install("msa") 
library(msa)
library(Biostrings) 

# Check current working directory
setwd("//wsl$/Ubuntu/home/rashedul/project/ProSEC/data/phylogeny")

# ---- ClustalW + UPGMA ----
seqs <- readDNAStringSet("./FFP/10-BetaSet_edited.fasta")
seqs <- readDNAStringSet("./FFP/17-BetaSet_edited.fasta")  # Use correct file
alignment <- msa(seqs, method = "ClustalW")
alignment_phangorn <- msaConvert(alignment, type = "phangorn::phyDat")
dist_matrix_clustalw <- dist.ml(alignment_phangorn)
tree_clustalw_upgma <- upgma(dist_matrix_clustalw)

# ---- ProSEQ + UPGMA ----
generate_upgma_tree <- function(matrix_csv, method_name) {
  corr_mat <- as.matrix(read.csv(matrix_csv, row.names = 1))
  dist_mat <- as.dist(1 - corr_mat)
  tree <- upgma(dist_mat)
  assign(paste0("tree_", method_name, "_upgma"), tree, envir = .GlobalEnv)
  return(tree)
}

# tree_proseq_upgma <- generate_upgma_tree("./FFP/10-BetaSet.proseq_matrix.csv", "proseq")
# test <- generate_upgma_tree("./FFP/10-BetaSet.esm2_small_matrix.csv", "esm2_small")
# pdf("../plots/esm2_small_UPGMA.pdf", width = 8, height = 6)
# plot(test, main = "esm2_small + UPGMA", cex = 0.8)
# dev.off()

# ---- Plotting the trees ----
pdf("../plots/clustalw_proseq_PLMs_UPGMA_Trees.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
plot(tree_clustalw_upgma, main = "ClustalW + UPGMA", cex = 0.8)
plot(tree_proseq_upgma, main = "ProSEQ + UPGMA", cex = 0.8)
dev.off()

matrix_files <- list.files("./FFP", pattern = "_matrix.csv$", full.names = TRUE)
method_names <- gsub("10-BetaSet.|\\.csv", "", basename(matrix_files))

pdf("../plots/all_methods_BetaSet_UPGMA_Trees.pdf", width = 16, height = 8)
par(mfrow = c(2, 3)) # Adjust rows/cols as needed

for (i in seq_along(matrix_files)) {
  tree <- generate_upgma_tree(matrix_files[i], method_names[i])
  plot(tree, main = paste(method_names[i], "+ UPGMA"), cex = 0.8)
}
dev.off()



