#!bin/bash

cd /home/rashedul/project/ProSEC/

# build db pkl 
for f in ./analysis/uniprot_data/uniprot_sprot*fasta; do
    echo "$f"
    filename=$(basename "$f")  
    python3 protsec.py \
        --fasta_path "$f" \
        --out_file "${filename}_db.pkl"
done 

# Annotate
for db in uni*_db.pkl; do
    echo "$db"
    db_name=$(basename "$db")  # Extract only the filename
    python3 annotate.py \
        --input_faa ./analysis/uniprot_data/uniprot_sprot_5000.fasta \
        --db "$db" \
        --dim 1024 \
        --out ./analysis/similarity_search/ProtSEC/"${db_name}_protsec_result.tsv"
done

rm *pkl

# Analysis 
cd ./analysis/similarity_search/ProtSEC/

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
done | paste - - | awk '{print $1 "\t" $2 "\t" 100*($2/5000)}' | sort -n -k3 >acc_percent_sigprot.txt
