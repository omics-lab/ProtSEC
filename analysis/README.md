## Analysis Directories

### similarity_search
Protein sequence similarity search benchmarking framework. Contains implementations and comparisons of different approaches including traditional BLAST, AI-embedding based methods, and ProtSEC. Includes benchmark scripts, dataset preparation tools, and evaluation pipelines for assessing search accuracy and performance across different protein databases.

**Key Scripts:**
- `sim_search_benchmark.sh` - Main benchmarking pipeline
- `dataset_preparation.py` - Dataset preprocessing
- `data_download.sh` - Download benchmark datasets
- `ncbi_nr_download.sh` - NCBI NR database download

```
bash ./similarity_search/sim_search_benchmark.sh
```

### EC_CARE_task1
Enzyme Commission (EC) classification analysis using k-nearest neighbors (k-NN) algorithm. Contains scripts for protein functional annotation based on precomputed distance matrices from various protein language models (PLMs). Includes training/test datasets, k-NN implementation in R, and batch evaluation workflows for comparing different embedding methods (ESM2, ProtBERT, ProtT5, ProtSEC).

**Key Scripts:**
- `knn.R` - Main k-NN classification script
- `knn_batch_eval.R` - Batch evaluation across methods
- `barplot_knn_accuracy.R` - Accuracy visualization

### phylogeny
Phylogenetic analysis comparing different protein similarity methods. Implements multiple sequence alignment (ClustalW) and UPGMA tree construction using various distance matrices from protein language models. Includes Branch Score Distance (BSD) analysis to evaluate how well different embedding methods recapitulate evolutionary relationships compared to traditional sequence alignment methods.

**Key Scripts:**
- `FFP/phylogeny.R` - Main phylogenetic analysis script with ClustalW, UPGMA trees, and BSD calculations

### uniprot_data
UniProt database processing and sequence extraction utilities. Contains scripts for downloading, filtering, and sampling protein sequences from UniProt SwissProt database in various dataset sizes (5K, 10K, 20K, 40K, 80K sequences). Includes tools for extracting specific protein families and preparing datasets for downstream analysis tasks.

**Key Scripts:**
- `extract_seq.py` - Sequence extraction and filtering utility

### Others

- Protein sequence similarity search using ProtSEC database file

Output `result.tsv` contains score in the 3rd column which is correlation value between query and hit.

```
cd ../
python3 annotate.py --input_faa ./data/QUERY.fasta \
    --db ./DB/mds_sms_db.pkl \
    --dim_reduct MDS \
    --dist_func SMS \
    --dim 1024 \
    --top_hit 1 \
    --out ./data/result.tsv
```


- Generate PLM Based Distance matrix 

```
cd ../
pip install biopython transformers torch sentencepiece
python3 get_plm_dist_mat_optimized.py -i data/phylogeny/FFP/17-BetaSet_edited.fasta -m esm2_small
```

- PLM Embedding
Code to run 'esm2_small', 'esm2_large', 'prot_bert', 'prot_t5' is available [here](https://github.com/Rajan-sust/GeneAnnotation) 


### Citation
Raju RS and Rashedul I. [ProtSEC: Ultrafast Protein Sequence Embedding in Complex Space Using Fast Fourier Transform. (2025)](https://www.biorxiv.org/content/10.1101/2025.08.17.670693v1). 