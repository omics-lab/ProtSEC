#!/bin/bash

## sample prep
cd /Users/rashedulislam/Documents/git_repos/ProtSEC/analysis/EC_CARE_task1/fastas/

# remove duplicated entry
seqkit rmdup backup/train_swissprot.fasta -o train_swissprot_unique.fasta

# downsampled testing data
cd ./analysis/EC_CARE_task1/fastas/ 
seqkit sample -n 5000 train_swissprot_unique.fasta -o sample_5000.fasta
seqkit sample -n 2000 train_swissprot_unique.fasta -o sample_2000.fasta
# prepare fasta for knn
cat sample_5000.fasta 30_protein_test.fasta 30-50_protein_test.fasta > sample_30-50_5000.fasta
cat sample_2000.fasta 30_protein_test.fasta 30-50_protein_test.fasta > sample_30-50_2000.fasta

# whole data (did not use)
cat train_swissprot_unique.fasta 30_protein_test.fasta 30-50_protein_test.fasta > train_test_swissprot_unique.fasta

## final 
#5000
python3 get_phase_dist_mat_optimized.py -n 1024 -i ./analysis/EC_CARE_task1/fastas/sample_30-50_5000.fasta -o ./analysis/EC_CARE_task1/fastas/sample_30-50_5000_ProtSEC.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_30-50_5000.fasta -m esm2_small -o ./analysis/EC_CARE_task1/fastas/sample_30-50_5000_esm2small.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_30-50_5000.fasta -m esm2_large -o ./analysis/EC_CARE_task1/fastas/sample_30-50_5000_esm2large.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_30-50_5000.fasta -m prot_bert -o ./analysis/EC_CARE_task1/fastas/sample_30-50_5000_protbert.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_30-50_5000.fasta -m prot_t5 -o ./analysis/EC_CARE_task1/fastas/sample_30-50_5000_prott5.csv

#2000
python3 get_phase_dist_mat_optimized.py -n 1024 -i ./analysis/EC_CARE_task1/fastas/sample_30-50_2000.fasta -o ./analysis/EC_CARE_task1/fastas/sample_30-50_2000_ProtSEC.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_30-50_2000.fasta -m esm2_small -o ./analysis/EC_CARE_task1/fastas/sample_30-50_2000_esm2small.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_30-50_2000.fasta -m esm2_large -o ./analysis/EC_CARE_task1/fastas/sample_30-50_2000_esm2large.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_30-50_2000.fasta -m prot_bert -o ./analysis/EC_CARE_task1/fastas/sample_30-50_2000_protbert.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_30-50_2000.fasta -m prot_t5 -o ./analysis/EC_CARE_task1/fastas/sample_30-50_2000_prott5.csv

#### backup codes
## ProtSEC
cd /Users/rashedulislam/Documents/git_repos/ProtSEC/
python3 get_phase_dist_mat.py -n 1024 -i ./analysis/EC_CARE_task1/fastas/sample_2000.fasta -o ./analysis/EC_CARE_task1/fastas/sample_2000_ProtSEC_matrix.csv
python3 get_plm_dist_mat.py -i ./analysis/EC_CARE_task1/fastas/sample_2000.fasta -m esm2_small -o ./analysis/EC_CARE_task1/fastas/sample_2000_esm2small_matrix.csv

## blastp
cd ./analysis/EC_CARE_task1/fastas/
makeblastdb -in sample_5000.fasta -dbtype prot -out sample_5000_db
blastp -query 30-50_protein_test.fasta -db sample_5000_db -outfmt 6 -out blastp_30-50.txt -max_target_seqs 1
