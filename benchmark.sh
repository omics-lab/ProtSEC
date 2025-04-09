#!bin/bash

## 5 datasets: 5k, 10k, 20k, 40k, 80k

## 1. blastp with 5 datasets 
# bash benchmark/blastp/blastp_analysis.sh

## 2. PA-SigPro benchmarking
# clone the repo from here (https://github.com/Rajan-sust/PA-SigPro-Pipeline)

# benchmark default sigprot with 5 datasets
# bash benchmark/PA-SigPro/PA-SigPro_analysis.sh 

# benchmark mds_tsne_umap dimension reduction menthods
# bash benchmark/mds_tsne_umap/mds_tsne_umap_analysis.sh

# benchmark `SMS`, `ASMP`, `SNN` distance methods
bash benchmark/SMS_ASMP_SNN/SMS_ASMP_SNN_analysis.sh 

## 3. AI-Embedding (prot_bert prot_t5 esm2)
# Clone the repo from here (https://github.com/Rajan-sust/GeneAnnotation) and run the benchmarking. 
# bash /mnt/c/GeneAnnotation/AI-Embedding_analysis.sh 

## 4. Performance benchmark 
# sigprot with 10k embedding performance
# /usr/bin/time -v benchmark/PA-SigPro/performance.sh &> benchmark/PA-SigPro/performance.log 

# AI methods 10k embedding performance
# /usr/bin/time -v benchmark/AI-Embedding/performance.sh &> benchmark/AI-Embedding/performance.log 
