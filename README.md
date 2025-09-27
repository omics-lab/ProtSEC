## ProtSEC architecture for protein seqence embedding
**ProtSEC** (**Prot**ein **S**equence **E**mbedding in **C**omplex Space) is an ultrafast method for embedding protein sequences using the discrete Fourier transform. Unlike large protein language models (PLMs), ProtSEC requires no training on sequence data. It is 20,000× faster and uses 85× less memory compared to the popular models like esm2_3B, esm2_35M, prot_t5 and prot_bert. ProtSEC is lightweight enough to run on personal or laptop computers, even for processing large protein sequence datasets. 

<p align="center">
  <img src="./analysis/figures_ProtSEC/github.jpeg" width="800"/>
</p>

### 1. Requirement

 - Python >= 3.10
 - Linux
 - macOS >= 13.5

### 2. Installation

 - Clone the repository and navigate to the project directory

```sh
git clone https://github.com/omics-lab/ProtSEC/
cd ProtSEC/
```

- Create a virtual environment and activate

```
python3 -m venv venv
source venv/bin/activate
```

- Upgrade `pip` and install the required dependencies

```
pip install --upgrade pip
pip install -r requirements.txt
```

### 3. Run ProtSEC

- Generate complex embedding using a FASTA file

```
# Available dimensionality reduction methods: `MDS`, `t-SNE`, `UMAP`
# Dist functions: `SMS`, `ASMP`, `SNN`

python3 protsec.py \
    --fasta_path ./data/DB.fasta \
    --dim_reduct MDS \
    --dist_func SMS \
    --dim 1024 \
    --db_dir_path ./DB \
    --db_filename mds_sms_db.pkl
```

- Protein sequence similarity search

Output `result.tsv` contains score in the 3rd column which is correlation value between query and hit.

```
python3 annotate.py --input_faa ./data/QUERY.fasta \
    --db ./DB/mds_sms_db.pkl \
    --dim_reduct MDS \
    --dist_func SMS \
    --dim 1024 \
    --top_hit 1 \
    --out ./data/result.tsv
```

- Generate phase correlation matrix using ProtSEC

`-n` : Dimension of the embedding. If you're working with a multi-FASTA file containing sequences of varying lengths, use the 75th percentile of sequence lengths. Otherwise, use the actual sequence length. Default is 1024.

```
python3 get_phase_dist_mat_optimized.py -n 1024 -i phosphatase.fa -o ProtSEC_matrix.csv
```

### 4. Benchmark 

- Benchmarking used in the manuscript

```
bash ./benchmark/benchmark.sh
```

- PLM Embedding
Code to run 'esm2_small', 'esm2_large', 'prot_bert', 'prot_t5' is available [here](https://github.com/Rajan-sust/GeneAnnotation) 

- Generate PLM Based Distance matrix

```
pip install biopython transformers torch sentencepiece
python3 get_plm_dist_mat_optimized.py -i data/phylogeny/FFP/17-BetaSet_edited.fasta -m esm2_small
```

### 5. Contact
Rashedul Islam, PhD (rashedul.gen@gmail.com)

### 6. Citation
Raju RS and Rashedul I. [ProtSEC: Ultrafast Protein Sequence Embedding in Complex Space Using Fast Fourier Transform. (2025)](https://www.biorxiv.org/content/10.1101/2025.08.17.670693v1).

### 7. License
Shield: [![CC BY-NC 4.0][cc-by-nc-shield]][cc-by-nc]

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg
