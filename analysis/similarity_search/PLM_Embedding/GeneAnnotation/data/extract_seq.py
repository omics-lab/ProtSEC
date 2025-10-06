import argparse
from Bio import SeqIO

def extract_sequences(input_fasta, output_fasta, search_string):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        # Parse the input FASTA file
        for record in SeqIO.parse(infile, "fasta"):
            # Check if the search string is in the sequence description (record.description)
            if search_string in record.description:
                # Write the matching sequence to the output file
                SeqIO.write(record, outfile, "fasta")
        print(f"Sequences matching '{search_string}' have been written to {output_fasta}")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Extract sequences from a FASTA file by matching the string in the description.")
    parser.add_argument('input_fasta', type=str, help="Input FASTA file")
    parser.add_argument('output_fasta', type=str, help="Output FASTA file to save matching sequences")
    parser.add_argument('search_string', type=str, help="String to match in the sequence descriptions")
    
    # Parse command-line arguments
    args = parser.parse_args()
    
    # Call the extract function with arguments
    extract_sequences(args.input_fasta, args.output_fasta, args.search_string)

if __name__ == "__main__":
    main()
