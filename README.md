## PA-SigPro-Pipeline
A bioinformatics pipeline for annotating protein sequences

### Installation

Clone the repository and navigate to the project directory:

```sh
git clone https://github.com/Rajan-sust/PA-SigPro-Pipeline
cd PA-SigPro-Pipeline/
```

Create a virtual environment and activate it:

```
python3 -m venv venv
source venv/bin/activate
```
Upgrade `pip` and install the required dependencies:
```
pip install --upgrade pip
pip install -r requirements.txt
```




###### Distance functions:

- Subtraction from Max Similarity (SMS)
```
D[i][j] = max(S) - S[i][j]
```

- Average of Self Minus Pairwise (ASMP)
```
D[i][j] = (S[i][i] + S[j][j]) / 2 - S[i][j]
```

- Shift and Normalize (SNN)
```
S' = -min(S) + S
D[i][j] = 1 - (S'[i][j] / max(S'))
```

### Usage
#### Building the Database

To build the database from a FASTA file, run:
```
# Available dimensionality reduction methods: `MDS`, `t-SNE`, `UMAP`
# Dist functions: `SMS`, `ASMP`, `SNN`

python3 db_build.py \
    --fasta_path ./data/DB.fasta \
    --dim_reduct MDS \
    --dist_func SMS \
    --db mds_sms_db.pkl
```

#### Annotating Sequences:

```
python3 annotate.py --input_faa ./data/QUERY.fasta \
    --db ./DB/mds_sms_db.pkl \
    --dim_reduct MDS \
    --dist_func SMS \
    --top_hit 1 \
    --out ./data/result.tsv
```

Output tsv contains score in the 3rd column which is correlation value between query and hit.

### Example data
Example data is provided in the `example_data` directory:


- `DB.fasta`: Example database FASTA file.
- `QUERY.fasta`: Example query FASTA file.


### Benchmark

Dependencies:
- blast
- SigProt

```
bash benchmark.sh
```

Code to run ProtBERT (Rostlab/prot_bert) and ESM2 (facebook/esm2) is available [here](https://github.com/Rajan-sust/GeneAnnotation) 

### Contact
Rashedul Islam, PhD (rashedul.gen@gmail.com)

### Citation
Raju RS and Rashedul I. SigProt: Ultra-fast protein sequence embedding method using evolutionary conservation
, (2025).

### License
This project is licensed under the MIT license.
