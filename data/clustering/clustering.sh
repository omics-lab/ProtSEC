#! /bin/bash

# clustering
cd /home/rashedul/project/ProSEC
pip install --upgrade pip
pip install -r requirements.txt

python3 get_phase_dist_mat.py -i data/clustering/protein_kinase/kinase.fa -o  data/clustering/protein_kinase/kinase.proseq_matrix.csv
python3 get_phase_dist_mat.py -i data/clustering/phosphatase/phosphatase.fa -o  data/clustering/phosphatase/phosphatase.proseq_matrix.csv
python3 get_phase_dist_mat.py -i data/clustering/protein_kinase/kinase.fa -o  data/clustering/protein_kinase/kinase.proseq_matrix.csv
python3 get_phase_dist_mat.py -i data/clustering/radical_sam/radicalsam.fa -o  data/clustering/radical_sam/radicalsam.proseq_matrix.csv
python3 get_phase_dist_mat.py -i data/clustering/deeploc/deeploc_data.fasta -o  data/clustering/deeploc/deeploc_data.proseq_matrix.csv

# dist mat plms
cd /home/rashedul/project/ProSEC

python3 get_plm_dist_mat.py -i data/clustering/deeploc/deeploc_data.fasta -m esm2_small -o data/clustering/deeploc/deeploc.esm2_small_matrix.csv
python3 get_plm_dist_mat.py -i data/clustering/deeploc/deeploc_data.fasta -m esm2_large -o data/clustering/deeploc/deeploc.esm2_large_matrix.csv
