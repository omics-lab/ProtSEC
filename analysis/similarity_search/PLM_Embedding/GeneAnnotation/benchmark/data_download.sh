#!bin/bash

cd /mnt/c/GeneAnnotation/data/

## Downlaod CUDASW benchmark datasets: https://github.com/asbschmidt/CUDASW4/tree/main

# small_db (uniprot sprot) (88 megabyte)
wget https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# medium_db (uniref50) (12 gigabyte)
wget https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta.gz

# large_db (uniprot trembl) (57 gigabyte)
wget -c https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

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

# spilt the query.fasta into three groups based on sizes of sequence
seqkit sort -l query.fasta | seqkit split -p 3 -O ../../../data/

cd /mnt/c/GeneAnnotation/data/

# db
small_db=uniprot_sprot.fasta
medium_db=uniref50.fasta
large_db=uniprot_trembl.fasta

# query
query6370=query.fasta
short_seq=stdin.part_001.fasta
medium_seq=stdin.part_002.fasta
long_seq=stdin.part_003.fasta
