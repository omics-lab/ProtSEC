#!bin/bash

# Analysis 
cd ./benchmark/AI-Embedding/

for f in *.tsv; do 
    echo "Processing: $f";
    awk '!seen[$1]++' "$f" | grep -v "Query" > "${f}.tophit";
done

k=uniprot_sprot_5k_against_10k_esm2.tsv.tophit
# count correct prediction
for k in *tophit; do
echo $k;
grep -oi 'GN=[^ ]*' $k | paste - - | awk '{if (tolower($1) == tolower($2)) print $0}' | wc -l; 
done | paste - - | awk '{print $1 "\t" $2 "\t" 100*($2/5000)}' | sort -n -k3 >acc_percent_AI_embd.txt
