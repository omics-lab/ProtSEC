#!/bin/bash

# Clone the repo from here (https://github.com/Rajan-sust/GeneAnnotation) and run the benchmarking. 

#cd /home/rashedul/project/GeneAnnotation/

cd /Users/rashedulislam/Documents/git_repos/GeneAnnotation

# Build database and run annotation
# for model in prot_bert prot_t5 esm2; do
for model in prot_t5; do  
    for f in data/uniprot_sprot_5000.fasta; do 
        echo "Processing file: $f with model: $model"
        filename=$(basename "$f" .fasta)  

        # Build database
        python3 db_build.py \
            --fasta_path "$f" \
            --collection "${filename}_${model}_tbl" \
            --model_name "$model" \
            --batch_size 2

        # Annotate with the same model
        python3 annotate.py \
            --input_faa "data/uniprot_sprot_5000.fasta" \
            --collection "${filename}_${model}_tbl" \
            --output_file ../PA-SigPro-Pipeline/benchmark/AI-Embedding/"${filename}_${model}_results.tsv" \
            --model_name "$model"
    done 
done


# # Acc Analysis 
# cd ../PA-SigPro-Pipeline/benchmark/AI-Embedding/

# for f in *.tsv; do 
#     echo "Processing: $f";
#     awk '!seen[$1]++' "$f" | grep -v "Query" > "${f}.tophit";
# done

# # k=uniprot_sprot_10000_esm2_results.tsv.tophit
# # count accurate prediction
# for k in *tophit; do
# echo $k;
# grep -oi 'GN=[^ ]*' $k | paste - - | awk '{if (tolower($1) == tolower($2)) print $0}' | wc -l; 
# done | paste - - | awk '{print $1 "\t" $2 "\t" 100*($2/5000)}' | sort -n -k3 >acc_percent_AI_embd.txt 