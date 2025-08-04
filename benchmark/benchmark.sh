#!bin/bash

cd /home/rashedul/project/ProSEC

## 1. blastp with 5 datasets 
bash benchmark/blastp/blastp_analysis.sh

## 2. PA-SigPro parameter benchmarking
# clone the repo from here (https://github.com/Rajan-sust/PA-SigPro-Pipeline) and install dependencies. 
# benchmark ProSEC with all parameters (mds_tsne_umap vs SMS_ASMP_SNN) for 5 datasets
# a DB directory is required
bash /home/rashedul/project/ProSEC/benchmark/dimreduct_distfunc/dimreduct_distfunc_analysis.sh 

## plot dimreduct_distfunc
# ERROR for packages if running from terminal
Rscript benchmark/dimreduct_distfunc/dimreduct_distfunc_plot.R

## 3. AI-Embedding (prot_bert prot_t5 esm2)
# Clone the repo from here (https://github.com/Rajan-sust/GeneAnnotation) and install dependencies. 
cd /mnt/c/GeneAnnotation/
bash /mnt/c/GeneAnnotation/benchmark/AI-Embedding_analysis.sh 

## 4. Computing performance benchmark [done]
# ProSEC with 5k embedding performance
cd /home/rashedul/project/ProSEC
/usr/bin/time -v benchmark/PA-SigPro/performance_evoprot.sh &> benchmark/PA-SigPro/performance_evoprot.log 

# PLM 5k embedding performance 
# run Docker and install python dependencies within GeneAnnotation
cd /mnt/c/GeneAnnotation/ 
bash benchmark/AI-Embedding/performance_esm2_protbert_prott5.sh  