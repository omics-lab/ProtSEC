import numpy as np



def cut_protein_sequence_center(protein_sequence, n):
    """
    Cut protein sequence from the center to achieve desired length n.
    
    If the sequence length is > n, it will remove the central portion
    and keep the N-terminal and C-terminal ends.
    
    Args:
        protein_sequence (str): Input protein sequence
        n (int): Desired output length
    
    Returns:
        str: Protein sequence of length n (if original length > n)
              or original sequence (if original length <= n)
    
    Examples:
        >>> cut_protein_sequence_center("ABCDEFGHIJ", 6)
        'ABCHIJ'
        >>> cut_protein_sequence_center("ABCDEFGHIJ", 4)
        'ABHIJ'
        >>> cut_protein_sequence_center("ABC", 6)
        'ABC'
    """
    # If sequence is already shorter than or equal to n, return as is
    if len(protein_sequence) <= n:
        return protein_sequence
    
    # Calculate how many residues to keep from each end
    total_to_keep = n
    
    # Split the kept residues equally between N-terminal and C-terminal
    # If odd number to keep, keep one more from the N-terminal
    n_terminal_keep = (total_to_keep) // 2
    c_terminal_keep = total_to_keep - n_terminal_keep
    
    # Extract N-terminal and C-terminal portions
    n_terminal = protein_sequence[:n_terminal_keep]
    c_terminal = protein_sequence[-c_terminal_keep:] if c_terminal_keep > 0 else ""
    
    # Combine N-terminal and C-terminal portions
    return n_terminal + c_terminal

def load_data_from_npz(file_path):
    """
    Load data from an NPZ file.
    
    Args:
        file_path (str): NPZ file to load
        
    Returns:
        tuple: (X_train, Y_train, X_test, Y_test)
    """
    try:
        data = np.load(file_path, allow_pickle=True)
        X_train = data['X_train']
        Y_train = data['Y_train']
        X_test = data['X_test']
        Y_test = data['Y_test']
        
        print(f"Data successfully loaded from '{file_path}'")
        return X_train, Y_train, X_test, Y_test
        
    except Exception as e:
        print(f"Error loading data from NPZ file: {e}")
        return None, None, None, None




if __name__ == "__main__":
    clf_data_path = 'Deeploc/clf_data.npz'
    X_train, Y_train, X_test, Y_test = load_data_from_npz(clf_data_path)
    if X_train is not None and Y_train is not None:
        print(f"Training data shape: {X_train.shape}, Labels shape: {Y_train.shape}")
        print(f"Test data shape: {X_test.shape}, Labels shape: {Y_test.shape}")
    else:
        print("Failed to load data from NPZ file.")

    # make protein embedder object
    