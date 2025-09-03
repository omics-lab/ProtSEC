import random
from Bio import SeqIO

# Set random seed for reproducibility
random.seed(42)

# File paths
input_fasta = "/mnt/c/GeneAnnotation/data/uniprot_sprot.fasta"
output_fasta = "/mnt/c/GeneAnnotation/data/uniprot_sprot_filtered.fasta"

# Words to filter out (case insensitive)
filter_words = {"putative", "uncharacter", "probable"}

# Function to check if a sequence should be kept
def should_keep(header):
    header_lower = header.lower()  # Convert header to lowercase for filtering
    contains_gn = "GN=" in header  # Case-sensitive check for "GN="
    excludes_words = not any(word in header_lower for word in filter_words)  # Exclude unwanted words
    return contains_gn and excludes_words  # Keep only if both conditions are met

# Read input, filter sequences, and save the filtered output
filtered_records = [
    record for record in SeqIO.parse(input_fasta, "fasta") if should_keep(record.description)
]

# Save the filtered sequences
with open(output_fasta, "w") as output_handle:
    SeqIO.write(filtered_records, output_handle, "fasta")

print(f"Filtering complete. Saved {len(filtered_records)} sequences to {output_fasta}")

# Sampling sizes
sample_sizes = [5000, 10000, 20000, 40000, 80000]

# Sample and write to different FASTA files
for sample_size in sample_sizes:
    if len(filtered_records) >= sample_size:  # Ensure enough sequences exist
        sampled_records = random.sample(filtered_records, sample_size)
        output_sample_file = f"/mnt/c/GeneAnnotation/data/uniprot_sprot_{sample_size}.fasta"
        with open(output_sample_file, "w") as sample_handle:
            SeqIO.write(sampled_records, sample_handle, "fasta")
        print(f"Saved {sample_size} sampled sequences to {output_sample_file}")
    else:
        print(f"Skipping {sample_size} sample (only {len(filtered_records)} sequences available)")

print("Sampling complete.")

