#!/bin/bash

# Clone the repo from here (https://github.com/Rajan-sust/GeneAnnotation) and run the benchmarking. 

GeneAnnotation=/mnt/c/GeneAnnotation
output=/home/rashedul/project/PA-SigPro-Pipeline/benchmark/AI-Embedding/

cd $GeneAnnotation

# for model in prot_bert prot_t5 esm2; do
for model in prot_bert prot_t5 ; do
    for f in data/uniprot_sprot_5000.fasta; do 
        echo "Processing file: $f with model: $model"
        filename=$(basename "$f" .fasta)  

        # Define log file name
        logfile="$output/${filename}_${model}_performance.txt"

        # Run and measure time/memory usage
        /usr/bin/time -v -o "$logfile" \
        python3 db_build.py \
            --fasta_path "$f" \
            --collection "${filename}_${model}_tbl" \
            --model_name "$model" \
            --batch_size 2

        echo "Log saved to $output"
    done 
done
