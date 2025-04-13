#!bin/bash

# /home/rashedul/project/PA-SigPro-Pipeline
# /mnt/c/GeneAnnotation/

## 1. blastp with 5 datasets 
# bash benchmark/blastp/blastp_analysis.sh

## 2. PA-SigPro parameter benchmarking
# clone the repo from here (https://github.com/Rajan-sust/PA-SigPro-Pipeline) and install dependencies. 
# benchmark sigprot with all parameters (mds_tsne_umap vs SMS_ASMP_SNN) for 5 datasets
# bash /home/rashedul/project/PA-SigPro-Pipeline/benchmark/dimreduct_distfunc/dimreduct_distfunc_analysis.sh 

# plot dimreduct_distfunc
# ERROR for packages
# Rscript benchmark/dimreduct_distfunc/dimreduct_distfunc_plot.R

## 3. AI-Embedding (prot_bert prot_t5 esm2)
# Clone the repo from here (https://github.com/Rajan-sust/GeneAnnotation) and install dependencies. 
# bash /mnt/c/GeneAnnotation/benchmark/AI-Embedding_analysis.sh 

## 4. Computing performance benchmark 
# sigprot with 10k embedding performance
# /usr/bin/time -v benchmark/PA-SigPro/performance.sh &> benchmark/PA-SigPro/performance.log 

# PLM 10k embedding performance
# /usr/bin/time -v benchmark/AI-Embedding/performance.sh &> benchmark/AI-Embedding/performance.log 
