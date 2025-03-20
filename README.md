### PA-SigPro-Pipeline
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

### Usage
###### Building the Database

To build the database from a FASTA file, run:
```
python3 db_build.py \
    --fasta_path ./data/DB.fasta \
    --db ./data/mmseq2_db.pkl
```

###### Annotating Sequences
```
python3 annotate.py \
    --input_faa ./data/QUERY.fasta \
    --db ./data/mmseq2_db.pkl \
    --out ./data/mmseq2_result.tsv
```


###### Example data
Example data is provided in the `example_data` directory:


- `DB.fasta`: Example database FASTA file.
- `QUERY.fasta`: Example query FASTA file

