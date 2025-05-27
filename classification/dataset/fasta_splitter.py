import sys
import random
from Bio import SeqIO
from io import StringIO

def split_fasta(input_stream, test_percent=0.2):
    """
    Split fasta file into test and training sets.
    
    Args:
        input_stream: File-like object containing FASTA data
        test_percent: Percentage of sequences to include in test set (0.05 = 5%)
        
    Returns:
        tuple: (test_sequences, train_sequences) as lists of SeqRecord objects
    """
    # Parse all sequences from the input
    sequences = list(SeqIO.parse(input_stream, "fasta"))
    
    # Determine number of sequences for test set
    total_sequences = len(sequences)
    test_count = max(1, int(total_sequences * test_percent))
    
    # Randomly select sequences for test set
    random.seed(31)  # Set seed for reproducibility
    test_indices = random.sample(range(total_sequences), test_count)
    test_indices_set = set(test_indices)
    
    # Split sequences into test and train sets
    test_sequences = [sequences[i] for i in test_indices]
    train_sequences = [sequences[i] for i in range(total_sequences) if i not in test_indices_set]
    
    return test_sequences, train_sequences

def main():
    # Read from stdin
    input_data = sys.stdin.read()
    
    # Split the data
    test_sequences, train_sequences = split_fasta(StringIO(input_data))
    
    # Write to output files
    SeqIO.write(test_sequences, "test.fasta", "fasta")
    SeqIO.write(train_sequences, "train.fasta", "fasta")
    
    # Print summary
    print(f"Split complete: {len(test_sequences)} sequences in test set, {len(train_sequences)} sequences in train set")

if __name__ == "__main__":
    main()
