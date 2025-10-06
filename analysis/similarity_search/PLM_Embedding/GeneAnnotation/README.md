### Prerequisite

- Prodigal (https://github.com/hyattpd/prodigal)
- HF Access Token (https://huggingface.co/docs/hub/en/security-tokens)

### Clone the repo
```
git clone https://github.com/Rajan-sust/GeneAnnotation.git
cd GeneAnnotation
```

### Python Environment Creation & Pkg Installation
```
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

### Build Vector DB of Protein

###### Features

- Support for multiple protein embedding models:
  - ProtBERT (Rostlab/prot_bert)
  - ESM2 (facebook/esm2)
  - ProtT5 (Rostlab)
  - OpenAI
- Multi-threaded processing of FASTA files
- Flexible and modular architecture

###### Command Line Arguments

- `--fasta_path`: Path to input FASTA file (required)
- `--collection`: Name of the collection to create in the Database (required)
- `--model_name`: Protein embedding model to use (choices: `prot_bert`, `esm2_small`, `esm2_large`, `prot_t5`, `openai`. default: `esm2`)
- `--batch_size`: Batch size for processing sequences (default: 50)



###### Usage
```
python3 db_build.py --fasta_path ./db.faa \
--collection esm2_tbl \
--model_name esm2_small \
--batch_size 2
```



### Gene Predictions
```
prodigal -i my.genome.fna  -g 11 -a protein.translations.faa
```

### Protein Annotation


###### Command Line Arguments

- `--input_faa`: Path to input FASTA file containing protein sequences (required)
- `--db_name`: Name of the Qdrant collection to search against (required)
- `--output_file`: Path to output TSV file for results (required)
- `--model_name`: Protein embedding model to use (choices: `prot_bert`, `esm2_small`, `esm2_large`, `prot_t5`, `openai`)


###### Example

```
python3 annotate.py --input_faa ./QUERY.fasta \
                    --collection esm2_tbl \
                    --output_file results_esm2.tsv \
                    --model_name esm2_small
```

###### Output Format

The tool generates a TSV file with the following columns:
- `Query_ID`: Identifier of the input sequence
- `Annotation`: Predicted protein annotation
- `Similarity_Score`: Similarity score `[-1.0, 1.0]` with the matched database entry
- `Status`: Processing status ('success', 'below_threshold', 'embedding_failed', or 'error')
