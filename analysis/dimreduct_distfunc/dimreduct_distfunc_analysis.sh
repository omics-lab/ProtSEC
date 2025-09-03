#!bin/bash

# build db pkl and run annotation
for dim in t-SNE MDS UMAP; do 
    for dis in SMS ASMP SNN; do
        echo "dim_reduct: $dim"
        echo "distance: $dis"
        
        for f in data/uniprot_sprot*fasta; do
            echo "Processing file: $f"
            filename=$(basename "$f" .fasta)

            # Build database
            python3 db_build.py \
                --fasta_path "$f" \
                --dim_reduct "$dim" \
                --dist_func "$dis" \
                --db_dir_path ./DB \
                --db_filename "${dim}_${dis}_${filename}_db.pkl"

            # Annotate 
            python3 annotate.py \
                --input_faa ./data/uniprot_sprot_5000.fasta \
                --db ./DB/"${dim}_${dis}_${filename}_db.pkl" \
                --dim_reduct "$dim" \
                --dist_func "$dis" \
                --top_hit 1 \
                --out ./benchmark/dimreduct_distfunc/"${dim}_${dis}_${filename}_analysis.tsv"

            # remove db    
            rm ./DB/"${dim}_${dis}_${filename}_db.pkl"
        done 
    done
done

# Analysis 
cd ./benchmark/dimreduct_distfunc/

# remove duplicated lines 
for f in *.tsv; do 
    echo "Processing: $f";
    awk '!seen[$1]++' "$f" > "${f}.tophit";
done

# count correct prediction
for k in *tophit; do
echo $k;
grep -oi 'GN=[^ ]*' $k | paste - - | awk '{if (tolower($1) == tolower($2)) print $0}' | wc -l; 
done | paste - - | awk '{print $1 "\t" $2 "\t" 100*($2/5000)}' >acc_percent_sigprot_methods.txt
