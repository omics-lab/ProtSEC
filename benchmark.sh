#!bin/bash

# # blastp 
# bash benchmark/blastp/blastp_analysis.sh

# # PA-SigPro 
# bash benchmark/PA-SigPro/PA-SigPro_analysis.sh 

# mds_tsne_umap
# bash benchmark/mds_tsne_umap/mds_tsne_umap_analysis.sh 

# AI-Embedding (need to add model generation)
# bash benchmark/blastp/AI-Embedding_analysis.sh 

## Performance 
# sigprot with 10k embedding
/usr/bin/time -v benchmark/PA-SigPro/performance.sh &> benchmark/PA-SigPro/performance.log 

# AI methods LLM
# /usr/bin/time -v benchmark/AI-Embedding/performance.sh &> benchmark/AI-Embedding/performance.log 
