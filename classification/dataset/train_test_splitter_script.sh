#!/bin/bash

# Clean up any existing output files
rm -f test_data.fasta train_data.fasta

# Process each UniRef fasta file
for file in U*.fasta; do
    echo "Processing $file..."
    
    # Sample 5% for test set &  95% for train set
    
    cat ${file} | python3 ../fasta_splitter.py
    # Append to the test file
    cat test.fasta >> test_data.fasta
    
    # Append to the train file
    cat train.fasta >> train_data.fasta
    
    # Clean up temporary files
    rm test.fasta train.fasta
done

# Report stats about the resulting files
echo "Done! Statistics:"
echo "Test set:"
seqkit stats test_data.fasta
echo "Training set:"
seqkit stats train_data.fasta
