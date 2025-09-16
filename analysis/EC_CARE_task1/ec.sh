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
python3 get_phase_dist_mat.py -n 1024 -i ./analysis/EC_CARE_task1/fastas/sample_2000.fasta -o ./analysis/EC_CARE_task1/fastas/sample_2000_ProtSEC_matrix.csv
python3 get_plm_dist_mat.py -i ./analysis/EC_CARE_task1/fastas/sample_2000.fasta -m esm2_small -o ./analysis/EC_CARE_task1/fastas/sample_2000_esm2small_matrix.csv
python3 get_plm_dist_mat.py -i ./analysis/EC_CARE_task1/fastas/sample_2000.fasta -m esm2_large -o ./analysis/EC_CARE_task1/fastas/sample_2000_esm2large_matrix.csv

# test performance
python3 get_phase_dist_mat_optimized.py -n 1024 -i ./analysis/EC_CARE_task1/fastas/sample_2000.fasta -o ./analysis/EC_CARE_task1/fastas/test.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_2000.fasta -m esm2_small -o ./analysis/EC_CARE_task1/fastas/test_esm2small_matrix.csv
