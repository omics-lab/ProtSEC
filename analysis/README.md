### Benchmark 

- Protein sequence similarity search

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

```
bash ./similarity_search/sim_search_benchmark.sh
```

- PLM Embedding
Code to run 'esm2_small', 'esm2_large', 'prot_bert', 'prot_t5' is available [here](https://github.com/Rajan-sust/GeneAnnotation) 

- Generate PLM Based Distance matrix

```
pip install biopython transformers torch sentencepiece
python3 get_plm_dist_mat_optimized.py -i data/phylogeny/FFP/17-BetaSet_edited.fasta -m esm2_small
```
