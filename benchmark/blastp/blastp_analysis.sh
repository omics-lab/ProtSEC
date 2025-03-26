#!bin/bash

cd ../data/

mkdir temp

# blastp db
for f in uniprot_sprot_[1-8]*.fasta; do 
    echo "Processing: $f"; 
    makeblastdb -in "$f" -dbtype prot -out temp/"${f%.fasta}_blastp_db";
done

# predict
for db in temp/uniprot_sprot_5000_blastp_db temp/uniprot_sprot_10000_blastp_db temp/uniprot_sprot_20000_blastp_db temp/uniprot_sprot_40000_blastp_db temp/uniprot_sprot_80000_blastp_db; do 
    echo "DB name: $db"

    # Extract database name without path for cleaner output filenames
    db_name=$(basename "$db")

    blastp -query "uniprot_sprot_5000.fasta" -db "$db" \
           -out "uniprot_sprot_5000_vs_${db_name}_blastp_results.txt" \
           -num_threads 4 \
           -outfmt "6 qseqid sseqid stitle evalue" \
           -max_target_seqs 1
done

# remove duplicated lines
for f in *_blastp_db_blastp_results.txt; do 
    echo $f; 
    awk '!seen[$1]++ && count[$1]++ == 0' $f >../benchmark/blastp/$f.duprem.txt;
done

# create id file of reference
less uniprot_sprot_5000.fasta | grep ">" | awk '{ gsub(">", ""); print $1 "\t" $0}' | awk '{ match($0, /GN=([^ ]+)/, arr); if (arr[1]) print arr[1] "\t" $0}' >../benchmark/blastp/uniprot_sprot_5000.fasta.id

rm *_blastp_db_blastp_results.txt
rm -r temp

# awk '{ match($0, /GN=([^ ]+)/, arr); if (arr[1]) print arr[1] "\t" $0}' uniprot_sprot_5000_vs_uniprot_sprot_10000_blastp_db_blastp_results.txt.duprem.txt | head
# k=uniprot_sprot_5000_vs_uniprot_sprot_10000_blastp_db_blastp_results.txt.duprem.txt

cd ../benchmark/blastp/

# count correct prediction
for k in *_db_blastp_results.txt.duprem.txt; do
echo $k;
comm -12 <(less uniprot_sprot_5000.fasta.id | awk '{print tolower($2)"|"tolower($1)}' | sort) <( less $k | awk '{ match($0, /GN=([^ ]+)/, arr); if (arr[1]) print arr[1] "\t" $0}' | awk '{print tolower($2)"|"tolower($1)}' | sort) | wc -l;
done | paste - - | awk '{print $1 "\t" $2 "\t" 100*($2/5000)}' | sort -k3 >acc_percent_blastp.txt
