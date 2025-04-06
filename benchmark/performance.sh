# provide full path
# cd benchmark
/usr/bin/time -v bash test.sh

/usr/bin/time -v python3 annotate.py ... 2> "./benchmark/mds_tsne_umap/time_${db_name}_${m}.log"

# Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.06
# Maximum resident set size (kbytes): 6244 is the peak RAM usage 
