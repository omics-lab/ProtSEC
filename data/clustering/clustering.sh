#! /bin/bash

# clustering
cd /home/rashedul/project/ProSEC
python3 gen-ph-cor.py -n 512 -i data/clustering/phosphatase/phosphatase.fa -o  data/clustering/phosphatase/phosphatase.proseq_matrix.csv
python3 gen-ph-cor.py -n 512 -i data/clustering/protein_kinase/kinase.fa -o  data/clustering/protein_kinase/kinase.proseq_matrix.csv
python3 gen-ph-cor.py -n 512 -i data/clustering/radical_sam/radicalsam.fa -o  data/clustering/radical_sam/radicalsam.proseq_matrix.csv
