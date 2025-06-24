from Bio import SeqIO
import numpy as np
from sklearn.utils import shuffle

subclass_to_parent = {
    'Cell.membrane-M': 'Cell membrane',
    'Cytoplasm-S': 'Cytoplasm',
    'Endoplasmic.reticulum-M': 'ER',
    'Endoplasmic.reticulum-S': 'ER',
    'Endoplasmic.reticulum-U': 'ER',
    'Extracellular-S': 'Extracellular',
    'Golgi.apparatus-M': 'Golgi apparatus',
    'Golgi.apparatus-S': 'Golgi apparatus',
    'Golgi.apparatus-U': 'Golgi apparatus',
    'Lysosome/Vacuole-M': 'Lysosome/Vacuole',
    'Lysosome/Vacuole-S': 'Lysosome/Vacuole',
    'Lysosome/Vacuole-U': 'Lysosome/Vacuole',
    'Mitochondrion-M': 'Mitochondrion',
    'Mitochondrion-S': 'Mitochondrion',
    'Mitochondrion-U': 'Mitochondrion',
    'Nucleus-M': 'Nucleus',
    'Nucleus-S': 'Nucleus',
    'Nucleus-U': 'Nucleus',
    'Peroxisome-M': 'Peroxisome',
    'Peroxisome-S': 'Peroxisome',
    'Peroxisome-U': 'Peroxisome',
    'Plastid-M': 'Plastid',
    'Plastid-S': 'Plastid',
    'Plastid-U': 'Plastid',
}


def load_and_process_data(fasta_file, random_state=42):
    """
    Load protein data from FASTA file and split into train/test sets.
    
    Args:
        fasta_file (str): Path to FASTA file
        random_state (int): Random seed for reproducibility
    
    Returns:
        tuple: (X_train, Y_train, X_test, Y_test) - shuffled arrays
    """
    X_train, Y_train, X_test, Y_test = [], [], [], []
    skipped_records = 0
    processed_records = 0
    
    subclass_keys = set(subclass_to_parent.keys())  # Use set for faster lookup
    
    print(f"Processing FASTA file: {fasta_file}")
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        description = record.description
        
        try:
            # Extract the subclass from the description
            # Assuming format like "Protein [Cell.membrane-M]" or similar
            parts = description.split()
            if len(parts) < 2:
                skipped_records += 1
                continue
                
            subclass = parts[1].strip()  # Remove brackets if present
            
            if subclass in subclass_keys:
                parent_class = subclass_to_parent[subclass]
                sequence = str(record.seq)
                
                # Check if this is a test record
                if 'test' in description.lower():
                    X_test.append(sequence)
                    Y_test.append(parent_class)
                else:
                    X_train.append(sequence)
                    Y_train.append(parent_class)
                
                processed_records += 1
            else:
                skipped_records += 1
                
        except (IndexError, AttributeError) as e:
            skipped_records += 1
            continue
    
    print(f"Processed {processed_records} records, skipped {skipped_records} records")
    print(f"Train set size: {len(X_train)}")
    print(f"Test set size: {len(X_test)}")
    
    # Convert to numpy arrays for efficient shuffling
    X_train = np.array(X_train, dtype=object)
    Y_train = np.array(Y_train, dtype=object)
    X_test = np.array(X_test, dtype=object)
    Y_test = np.array(Y_test, dtype=object)
    
    # Shuffle train and test sets separately while maintaining X-Y correspondence
    if len(X_train) > 0:
        X_train, Y_train = shuffle(X_train, Y_train, random_state=random_state)
        print("Training data shuffled")
    
    if len(X_test) > 0:
        X_test, Y_test = shuffle(X_test, Y_test, random_state=random_state)
        print("Test data shuffled")
    
    return X_train, Y_train, X_test, Y_test


def save_data_to_npz(X_train, Y_train, X_test, Y_test, filename='protein_data.npz'):
    """
    Save the processed data arrays to an NPZ file.
    
    Args:
        X_train, Y_train, X_test, Y_test: Data arrays to save
        filename (str): Output filename for the NPZ file
    """
    try:
        np.savez_compressed(
            filename,
            X_train=X_train,
            Y_train=Y_train,
            X_test=X_test,
            Y_test=Y_test
        )
        print(f"\nData successfully saved to '{filename}'")
        print(f"File size: {get_file_size(filename):.2f} MB")
        
    except Exception as e:
        print(f"Error saving data to NPZ file: {e}")


def load_data_from_npz(filename):
    """
    Load data from an NPZ file.
    
    Args:
        filename (str): NPZ file to load
        
    Returns:
        tuple: (X_train, Y_train, X_test, Y_test)
    """
    try:
        data = np.load(filename, allow_pickle=True)
        X_train = data['X_train']
        Y_train = data['Y_train']
        X_test = data['X_test']
        Y_test = data['Y_test']
        
        print(f"Data successfully loaded from '{filename}'")
        return X_train, Y_train, X_test, Y_test
        
    except Exception as e:
        print(f"Error loading data from NPZ file: {e}")
        return None, None, None, None


def get_file_size(filename):
    """Get file size in MB."""
    import os
    try:
        size_bytes = os.path.getsize(filename)
        return size_bytes / (1024 * 1024)  # Convert to MB
    except:
        return 0


def print_data_summary(X_train, Y_train, X_test, Y_test):
    """Print summary statistics of the loaded data."""
    from collections import Counter
    
    print("\n=== Data Summary ===")
    print(f"Training samples: {len(X_train)}")
    print(f"Test samples: {len(X_test)}")
    
    if len(Y_train) > 0:
        train_counts = Counter(Y_train)
        print(f"\nTraining set class distribution:")
        for class_name, count in sorted(train_counts.items()):
            print(f"  {class_name}: {count}")
    
    if len(Y_test) > 0:
        test_counts = Counter(Y_test)
        print(f"\nTest set class distribution:")
        for class_name, count in sorted(test_counts.items()):
            print(f"  {class_name}: {count}")
    
    if len(X_train) > 0:
        avg_train_length = np.mean([len(seq) for seq in X_train])
        print(f"\nAverage training sequence length: {avg_train_length:.2f}")
    
    if len(X_test) > 0:
        avg_test_length = np.mean([len(seq) for seq in X_test])
        print(f"Average test sequence length: {avg_test_length:.2f}")


if __name__ == '__main__':
    # Main execution
    fasta_file = 'Deeploc/deeploc_data.fasta'  # Replace with your actual FASTA file path
    npz_filename = 'Deeploc/clf_data.npz'  # Output NPZ filename
    
    try:
        # Load and process data with shuffling
        X_train, Y_train, X_test, Y_test = load_and_process_data(
            fasta_file, 
            random_state=31  # For reproducible results
        )
        
        # Print summary
        print_data_summary(X_train, Y_train, X_test, Y_test)
        
        # Save data to NPZ file
        save_data_to_npz(X_train, Y_train, X_test, Y_test, npz_filename)
        
        # Demonstrate loading the saved data
        print(f"\n=== Loading saved data ===")
        X_train_loaded, Y_train_loaded, X_test_loaded, Y_test_loaded = load_data_from_npz(npz_filename)
        
        if X_train_loaded is not None:
            print(f"Loaded training samples: {len(X_train_loaded)}")
            print(f"Loaded test samples: {len(X_test_loaded)}")
            
            # Verify data integrity
            if np.array_equal(X_train, X_train_loaded) and np.array_equal(Y_train, Y_train_loaded):
                print("✓ Data integrity verified - saved and loaded data match")
            else:
                print("✗ Warning: Loaded data differs from original")
            
    except FileNotFoundError:
        print(f"Error: Could not find FASTA file '{fasta_file}'")
        print("Please check the file path and try again.")
    except Exception as e:
        print(f"An error occurred: {e}")