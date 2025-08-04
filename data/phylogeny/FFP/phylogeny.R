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
corr_mat <- as.matrix(read.csv("FFP/10-BetaSet.proseq_matrix.csv", row.names = 1))
corr_mat <- as.matrix(read.csv("FFP/17-BetaSet.proseq_matrix.csv", row.names = 1))
dist_mat_proseq <- as.dist(1 - corr_mat)
tree_proseq_upgma <- upgma(dist_mat_proseq)

# ---- Save both trees to one PDF ----
pdf("../../plots/Clustalw_proseq_UPGMA_Trees.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
plot(tree_clustalw_upgma, main = "ClustalW + UPGMA", cex = 0.8)
plot(tree_proseq_upgma, main = "ProSEQ + UPGMA", cex = 0.8)
dev.off()

