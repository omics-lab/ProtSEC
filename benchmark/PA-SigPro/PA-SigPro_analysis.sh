#!bin/bash

# mkdir data/temp
cd ../

for f in data/uniprot_sprot*fasta; do
    echo "$f"
    filename=$(basename "$f")  # Extract only the filename
    python3 db_build.py \
        --fasta_path "$f" \
        --db "${filename}_db.pkl"
done

# Annotate
for db in DB/uni*_db.pkl; do
    echo "$db"
    db_name=$(basename "$db")  # Extract only the filename
    python3 annotate.py \
        --input_faa ./data/uniprot_sprot_5000.fasta \
        --db "$db" \
        --out ./benchmark/PA-SigPro/"${db_name}_PA-SigPro_result.tsv"
done

rm -r DB

# Analysis 
cd ./benchmark/PA-SigPro

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
