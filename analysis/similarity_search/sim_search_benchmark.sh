#!bin/bash

cd ProSEC/analysis/similarity_search/

## 1. blastp with 5 datasets 
bash blastp/blastp_analysis.sh

## 2. ProtSEC parameter benchmarking
# benchmark ProSEC with all parameters (mds_tsne_umap vs SMS_ASMP_SNN) for 5 datasets
bash dimreduct_distfunc/dimreduct_distfunc_analysis.sh 

## plot dimreduct_distfunc
# ERROR for packages if running from terminal
Rscript dimreduct_distfunc/dimreduct_distfunc_plot.R

## 3. AI-Embedding (prot_bert prot_t5 esm2)
# Clone the repo from here (https://github.com/Rajan-sust/GeneAnnotation) and install dependencies. 
bash GeneAnnotation/benchmark/plm.sh 

## 4. Computing performance benchmark [done]
# ProSEC with 5k embedding performance
/usr/bin/time -v ./ProSEC/performance_protsec.sh &> ./ProSEC/performance_protsec.log

# PLM 5k embedding performance 
# run Docker and install python dependencies within GeneAnnotation
bash GeneAnnotation/benchmark/AI-Embedding/performance_esm2_protbert_prott5.sh