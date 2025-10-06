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


# blastp
cd /Users/rashedulislam/Documents/git_repos/GeneAnnotation/data

# blastp with defaul parameters
makeblastdb -in DB.fasta -dbtype prot -out my_blast_db
# output 1 and 10 seqs
blastp -query QUERY.fasta -db my_blast_db -out blastp_results.txt -num_threads 4 -evalue 1e-10 -outfmt "6 qseqid sseqid stitle evalue" -max_target_seqs 1
blastp -query QUERY.fasta -db my_blast_db -out blastp_results_10.txt -num_threads 4 -evalue 1e-10 -outfmt "6 qseqid sseqid stitle evalue" -max_target_seqs 10

blastp -query ATP-dependent_DNA_helicase_RecG_only.fasta.gz -db my_blast_db -out blastp_results_ATP.txt -num_threads 4 -evalue 1e-10 -outfmt "6 qseqid sseqid stitle evalue" -max_target_seqs 1

awk '!seen[$1]++' blastp_results.txt >blastp_results_remdupli.txt


# downlaod from: https://www.uniprot.org/uniprotkb?query=ATP-dependent+DNA+helicase+RecG
uniprotkb_ATP_dependent_DNA_helicase_Re_2025_03_09.fasta.gz

# extract only "ATP-dependent DNA helicase RecG"
python extract_seq.py uniprotkb_ATP_dependent_DNA_helicase_Re_2025_03_09.fasta ATP-dependent_DNA_helicase_RecG_only.fasta "ATP-dependent DNA helicase RecG OS="


cd /Users/rashedulislam/Documents/git_repos/GeneAnnotation

## mmseq2
# single line code
mmseqs easy-search ATP-dependent_DNA_helicase_RecG_only.fasta DB.fasta alnRes.m8 tmp

# prebuild db
mmseqs createdb examples/DB.fasta targetDB
mmseqs createindex targetDB tmptrfgeAF	
mmseqs easy-search examples/QUERY.fasta targetDB alnRes.m8 tmp