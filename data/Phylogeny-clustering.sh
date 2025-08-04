# data download
esearch -db protein -query "P02024" | efetch -format fasta  >>10-BetaSet.fasta
esearch -db protein -query "ALU64020" | efetch -format fasta  >>17-BetaSet.fasta

# phylogeny
cd /home/rashedul/project/ProSEC
python3 gen-ph-cor.py -n 512 -i data/phylogeny/FFP/10-BetaSet_edited.fasta -o  data/phylogeny/FFP/10-BetaSet.proseq_matrix.csv
python3 gen-ph-cor.py -n 512 -i data/phylogeny/FFP/17-BetaSet_edited.fasta -o  data/phylogeny/FFP/17-BetaSet.proseq_matrix.csv

# clustering
cd /home/rashedul/project/ProSEC
python3 gen-ph-cor.py -n 512 -i data/clustering/phosphatase/phosphatase.fa -o  data/clustering/phosphatase/phosphatase.proseq_matrix.csv
python3 gen-ph-cor.py -n 512 -i data/clustering/protein_kinase/kinase.fa -o  data/clustering/protein_kinase/kinase.proseq_matrix.csv
python3 gen-ph-cor.py -n 512 -i data/clustering/radical_sam/radicalsam.fa -o  data/clustering/radical_sam/radicalsam.proseq_matrix.csv