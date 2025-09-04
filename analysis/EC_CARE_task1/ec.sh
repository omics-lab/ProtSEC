#!/bin/bash

## sample prep
cd /Users/rashedulislam/Documents/git_repos/ProtSEC/analysis/EC_CARE_task1/fastas/

# remove duplicated entry
seqkit rmdup backup/train_swissprot.fasta -o train_swissprot_unique.fasta

# downsampled testing data
seqkit sample -n 2000 train_swissprot_unique.fasta -o temp.fasta
cat temp.fasta 30_protein_test.fasta 30-50_protein_test.fasta > sample_2000.fasta
rm temp.fasta

# whole data
cat train_swissprot_unique.fasta 30_protein_test.fasta 30-50_protein_test.fasta > train_test_swissprot_unique.fasta

## ProtSEC
cd /Users/rashedulislam/Documents/git_repos/ProtSEC/

# downsampled testing data
python3 get_phase_dist_mat.py -n 1024 -i ./analysis/EC_CARE_task1/fastas/sample_2000.fasta -o ./analysis/EC_CARE_task1/fastas/sample_2000_ProtSEC_matrix.csv

# whole data

