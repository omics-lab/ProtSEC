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
- Minimum of 10GB disk space is required

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
- Time: ~2 minute on a personal computer 

```
pip install --upgrade pip
pip install -r requirements.txt
```

### 3. Run ProtSEC

- Generate complex embedding using a FASTA file
- Time: ~1 second on a personal computer
- Output file: `protein_embedding_ProtSEC.pkl` (size ~80MB)

```
python3 protsec.py \
    --fasta_path ./analysis/uniprot_data/uniprot_sprot_5000.fasta
```

#### Usage

```
python protsec.py -h
usage: protsec.py [-h] --fasta_path FASTA_PATH [--dim DIM] [--num_threads NUM_THREADS] [--dim_reduct {UMAP,t-SNE,MDS}] [--dist_func {SMS,ASMP,SNN}]
                  [--out_file OUT_FILE]

Build protein vector database from FASTA file

options:
  -h, --help            show this help message and exit
  --fasta_path FASTA_PATH
                        Path to input FASTA file (default: None)
  --dim DIM             Dimensionality of the embeddings (default: 1024)
  --num_threads NUM_THREADS
                        Number of worker threads to use (default: 1/4 of CPU cores) (default: 6)
  --dim_reduct {UMAP,t-SNE,MDS}
                        Algorithm for dimensionality reduction (default: MDS)
  --dist_func {SMS,ASMP,SNN}
                        Distance function for computing distance (default: ASMP)
  --out_file OUT_FILE   Output file path for the embeddings (default: protein_embedding_ProtSEC.pkl)
```

#### Generate phase correlation matrix using ProtSEC

- To select the optimum dimension of the embedding `n`, please check the Methods section of the paper. 
- Default `-n` is 1024.

```
fasta=./analysis/uniprot_datauniprot_sprot_5000.fasta
python3 get_phase_dist_mat_optimized.py -n 1024 -i $fasta -o ProtSEC_matrix.csv
```

### 4. Contact
Rashedul Islam, PhD (rashedul.gen@gmail.com)

### 5. Citation
Raju RS and Rashedul I. [ProtSEC: Ultrafast Protein Sequence Embedding in Complex Space Using Fast Fourier Transform. (2025)](https://www.biorxiv.org/content/10.1101/2025.08.17.670693v1).

### 6. License
[![CC BY-NC 4.0][cc-by-nc-shield]][cc-by-nc]

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg
