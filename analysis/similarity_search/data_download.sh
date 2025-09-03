#!bin/bash

cd /mnt/c/GeneAnnotation/data/

## Downlaod CUDASW benchmark datasets: https://github.com/asbschmidt/CUDASW4/tree/main

# small_db (uniprot sprot) (88 megabyte)
wget https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# medium_db (uniref50) (12 gigabyte)
wget https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta.gz

# large_db (uniprot trembl) (57 gigabyte)
wget https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

# 20 query sequences used in CUDASW
wget https://github.com/asbschmidt/CUDASW4/allqueries.fasta

## Download Chorus benchmark datasets: https://github.com/Bio-Acc/Chorus-Reproduction-Scripts?tab=readme-ov-file#download-data

wget http://wwwuser.gwdg.de/~compbiol/mmseqs2/mmseqs2-benchmark.tar.gz
tar -zxvf mmseqs2-benchmark.tar.gz
rm -r mmseqs2-benchmark.tar.gz 

cd /mnt/c/GeneAnnotation/data/mmseqs2-benchmark-pub/db
# query.fasta: randoamly selected 6370 seqs from SCOP database (version 1.75)  
# query100x.fasta: 637000 seqs from query.fasta
# targetdb.fasta: UniRef50 (June 2015 version)

# spilt the query.fasta into three groups base on sizes of sequence

seqkit sort -l query.fasta | seqkit split -p 3 -O ../../../data/

cd /mnt/c/GeneAnnotation/data/

# db
small_db=uniprot_sprot.fasta
medium_db=uniref50.fasta
large_db=uniprot_trembl.fasta

# query
small_seq=stdin.part_001.fasta
medium_seq=stdin.part_002.fasta
long_seq=stdin.part_003.fasta

# uniref50 keep only with GN id
from Bio import SeqIO

# Define input and output file paths
input_fasta = "/mnt/c/GeneAnnotation/data/uniprot_sprot.fasta"  # Replace with your actual input file
output_fasta = "/mnt/c/GeneAnnotation/data/uniprot_sprot_filtered.fasta"

# Words to filter out (case insensitive)
filter_words = {"putative", "uncharacter", "probab"}  # Stored in lowercase

# Function to check if a sequence should be kept
def should_keep(header):
    header_lower = header.lower()  # Convert header to lowercase for filtering
    contains_gn = "GN=" in header  # Case-sensitive check for "GN="
    excludes_words = not any(word in header_lower for word in filter_words)  # Check filter words
    return contains_gn and excludes_words  # Keep only if it meets both conditions

# Read input and write filtered sequences
with open(output_fasta, "w") as output_handle:
    filtered_records = (
        record for record in SeqIO.parse(input_fasta, "fasta") if should_keep(record.description)
    )
    SeqIO.write(filtered_records, output_handle, "fasta")

print(f"Filtering complete. Saved filtered sequences to {output_fasta}")

