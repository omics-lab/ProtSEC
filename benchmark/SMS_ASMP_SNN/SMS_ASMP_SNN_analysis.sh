#!bin/bash

# build db pkl and run annotation
for m in SMS ASMP SNN; do 
    echo "Method: $m"
    for f in data/uniprot_sprot*fasta; do
        echo "Processing file: $f"
        filename=$(basename "$f")

        # build db 
        python3 db_build.py \
            --fasta_path "$f" \
            --dim_reduct MDS \
            --dist_func $m \
            --db "${m}_${filename}_db.pkl"

        # annotate    
        python3 annotate.py \
            --input_faa ./data/uniprot_sprot_5000.fasta \
            --db ./DB/"${m}_${filename}_db.pkl" \
            --dim_reduct MDS \
            --dist_func $m \
            --top_hit 1 \
            --out ./benchmark/SMS_ASMP_SNN/"${m}_${filename}_SigProt_result.tsv"
    done 
done

# rm -r DB

# Analysis 
cd ./benchmark/SMS_ASMP_SNN/

# remove duplicated lines 
for f in *.tsv; do 
    echo "Processing: $f";
    awk '!seen[$1]++' "$f" > "${f}.tophit";
done

# k=uniprot_sprot_10000.fasta_db.pkl_PA-SigPro_result.tsv.tophit
# count correct prediction
for k in *tophit; do
echo $k;
grep -oi 'GN=[^ ]*' $k | paste - - | awk '{if (tolower($1) == tolower($2)) print $0}' | wc -l; 
done | paste - - | awk '{print $1 "\t" $2 "\t" 100*($2/5000)}' >acc_percent_sigprot_methods.txt
