#!bin/bash
cd ../../
python3 get_phase_dist_mat_optimized.py -n 1024 -i ./analysis/EC_CARE_task1/fastas/sample_30-50_2000.fasta -o ./analysis/EC_CARE_task1/fastas/sample_30-50_2000_ProtSEC.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_30-50_2000.fasta -m esm2_small -o ./analysis/EC_CARE_task1/fastas/sample_30-50_2000_esm2small.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_30-50_2000.fasta -m esm2_large -o ./analysis/EC_CARE_task1/fastas/sample_30-50_2000_esm2large.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_30-50_2000.fasta -m prot_bert -o ./analysis/EC_CARE_task1/fastas/sample_30-50_2000_protbert.csv
python3 get_plm_dist_mat_optimized.py -i ./analysis/EC_CARE_task1/fastas/sample_30-50_2000.fasta -m prot_t5 -o ./analysis/EC_CARE_task1/fastas/sample_30-50_2000_prott5.csv

