cd /Users/rashedulislam/Documents/git_repos/PA-SigPro-Pipeline

# db
python3 annotate.py \
    --input_faa ./example_data/QUERY.fasta \                                       
    --db ./DB/mmseq2_db.pkl \
    --out mmseq2_result.tsv     

# predict
python3 annotate.py \
    --input_faa ../GeneAnnotation/data/ATP-dependent_DNA_helicase_RecG_only.fasta \
    --db ./DB/mmseq2_db.pkl \
    --out ATP-dependent_DNA_helicase_RecG_results.tsv

less ATP-dependent_DNA_helicase_RecG_results.tsv | awk '!seen[$1]++' | awk -F'\t' '{print $2}' | wc -l
# 34392
less ATP-dependent_DNA_helicase_RecG_results.tsv | awk '!seen[$1]++' | awk -F'\t' '{print $2}' | grep "GN=recG" | wc -l
# 25196

# 25196/34392 73% acc

# add corr cutoff >.2
less ATP-dependent_DNA_helicase_RecG_results.tsv | awk '!seen[$1]++' | awk -F'\t' '$3>.2{print $2}' | grep "ATP-dependent DNA helicase RecG OS" | wc -l 

cd /Users/rashedulislam/Documents/git_repos/GeneAnnotation/data

## mmseq2
# single line code
mmseqs easy-search ATP-dependent_DNA_helicase_RecG_only.fasta DB.fasta alnRes.m8 tmp

# # prebuild db
# mmseqs createdb examples/DB.fasta targetDB
# mmseqs createindex targetDB tmp
# mmseqs easy-search examples/QUERY.fasta targetDB alnRes.m8 tmp

less alnRes.m8 | awk '!seen[$1]++' | awk '{print $2}'  >mmseq_hits.txt
less DB.fasta| grep ">" >mmseq_db.txt 

while read line; do grep $line mmseq_db.txt; done <mmseq_hits.txt | grep "GN=recG" >mmseq_hitsname.txt 

wc -l mmseq_hits.txt
# 31581
wc -l mmseq_hitsname.txt
# 31549 

# total ATP query = 34392
# acc 31549/34392 = 91.73%

# blastp
cd /Users/rashedulislam/Documents/git_repos/GeneAnnotation/data

# blastp with defaul parameters
makeblastdb -in DB.fasta -dbtype prot -out my_blast_db
# output 1 seq
blastp -query ATP-dependent_DNA_helicase_RecG_only.fasta -db my_blast_db -out blastp_results.txt -num_threads 4 -outfmt "6 qseqid sseqid stitle evalue" -max_target_seqs 1

awk '!seen[$1]++' blastp_results.txt | wc -l     
# 34361
awk '!seen[$1]++' blastp_results.txt | grep "GN=recG" | wc -l 
# 31571

# total ATP query = 34392
# acc 31571/34392 = 91.80%
