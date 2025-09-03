#! /bin/bash

# data download
esearch -db protein -query "P02024" | efetch -format fasta  >>10-BetaSet.fasta
esearch -db protein -query "ALU64020" | efetch -format fasta  >>17-BetaSet.fasta

# phylogeny
cd /home/rashedul/project/ProSEC
python3 gen-ph-cor.py -n 512 -i data/phylogeny/FFP/10-BetaSet_edited.fasta -o  data/phylogeny/FFP/10-BetaSet.proseq_matrix.csv
python3 gen-ph-cor.py -n 512 -i data/phylogeny/FFP/17-BetaSet_edited.fasta -o  data/phylogeny/FFP/17-BetaSet.proseq_matrix.csv
