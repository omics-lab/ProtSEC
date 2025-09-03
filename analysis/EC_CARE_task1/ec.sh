# ec 
cd /Users/rashedulislam/Documents/git_repos/ProtSEC/

# Loop through all FASTA files in the fastas directory
for fasta_file in ./data/EC/task1/fastas/*.fasta; do
    # Extract filename without path and extension
    filename=$(basename "$fasta_file" .fasta)
    
    # Create output filename
    output_file="./data/EC/task1/fastas/${filename}_score_matrix.csv"
    
    echo "Processing: $fasta_file"
    echo "Output: $output_file"
    
    # Run get_phase_dist_mat.py
    python3 get_phase_dist_mat.py -n 1024 -i "$fasta_file" -o "$output_file"
    
    echo "Completed: $filename"
    echo "---"
done

echo "All files processed!"

#python3 get_phase_dist_mat.py -n 1024 -i ./data/EC/task1/fastas/train_swissprot.fasta -o ./data/EC/task1/fastas/train_swissprot_score_matrix.csv