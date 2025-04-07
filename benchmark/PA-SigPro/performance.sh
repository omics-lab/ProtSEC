# build db pkl 
for f in data/uniprot_sprot_10000.fasta; do
    echo "$f"
    filename=$(basename "$f")  
    python3 db_build.py \
        --fasta_path "$f" \
        --db "${filename}_db.pkl"
done 

# Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.06
# Maximum resident set size (kbytes): 6244 is the peak RAM usage 