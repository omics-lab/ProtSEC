from embedder import ProteinEmbedder
from Bio import SeqIO

import numpy as np
from scipy.fft import ifft
import argparse
import pandas as pd
from multiprocessing import Pool
import os
import time


def calculate_distance(F1, F2):
    """Optimized single distance calculation"""
    cross_power_spectrum = F1 * np.conj(F2)
    abs_cross_power = np.abs(cross_power_spectrum)
    # Avoid division by zero and use in-place operation
    np.divide(cross_power_spectrum, abs_cross_power, out=cross_power_spectrum, where=(abs_cross_power != 0))
    correlation = ifft(cross_power_spectrum)
    peak_value = np.max(np.abs(correlation))
    distance = 1 - peak_value
    return distance


def calculate_distance_vectorized_row(F1, F_all):
    """Vectorized calculation for one row against all sequences"""
    n = len(F_all)
    distances = np.zeros(n)
    
    # Broadcast F1 against all sequences
    F1_broadcast = np.tile(F1, (n, 1))  # Shape: (n, fft_dim)
    F_all_array = np.array(F_all)       # Shape: (n, fft_dim)
    
    # Vectorized cross-power spectrum calculation
    cross_power = F1_broadcast * np.conj(F_all_array)
    
    # Avoid division by zero
    abs_cross_power = np.abs(cross_power)
    abs_cross_power[abs_cross_power == 0] = 1e-12
    cross_power /= abs_cross_power
    
    # Vectorized IFFT
    correlations = ifft(cross_power, axis=1)
    
    # Find peak values for each correlation
    peak_values = np.max(np.abs(correlations), axis=1)
    
    # Convert to distances
    distances = 1 - peak_values
    
    return distances


def process_chunk(args):
    """Process a chunk of rows for multiprocessing"""
    i_start, i_end, features, n, chunk_id = args
    
    print(f"Processing chunk {chunk_id}: rows {i_start}-{i_end-1}")
    
    chunk_size = i_end - i_start
    chunk_distances = np.zeros((chunk_size, n))
    
    for local_i, i in enumerate(range(i_start, i_end)):
        if i == 0 or (i + 1) % 10 == 0:
            print(f"Chunk {chunk_id}: processing row {i+1}/{n}")
            
        # Use vectorized calculation for this row
        row_distances = calculate_distance_vectorized_row(features[i], features)
        chunk_distances[local_i] = np.round(row_distances, 5)
        
        # Set diagonal to 0 (distance to self)
        chunk_distances[local_i, i] = 0.0
    
    return i_start, i_end, chunk_distances


def calculate_distances_optimized(features, num_processes=None):
    """Optimized pairwise distance calculation using multiprocessing and vectorization"""
    n = len(features)
    
    if num_processes is None:
        num_processes = min(os.cpu_count(), 8)  # Limit to 8 processes to avoid memory issues
    
    print(f"Using {num_processes} processes for calculation")
    
    # Calculate optimal chunk size
    chunk_size = max(1, n // (num_processes * 2))
    print(f"Chunk size: {chunk_size}")
    
    # Create chunks
    chunks = []
    for i, start in enumerate(range(0, n, chunk_size)):
        end = min(start + chunk_size, n)
        chunks.append((start, end, features, n, i))
    
    print(f"Created {len(chunks)} chunks")
    
    # Process chunks in parallel
    start_time = time.time()
    
    if num_processes == 1:
        # Sequential processing for debugging
        results = [process_chunk(chunk) for chunk in chunks]
    else:
        # Parallel processing
        with Pool(processes=num_processes) as pool:
            results = pool.map(process_chunk, chunks)
    
    processing_time = time.time() - start_time
    print(f"Chunk processing completed in {processing_time:.2f} seconds")
    
    # Assemble results
    print("Assembling distance matrix...")
    distances = np.zeros((n, n), dtype=np.float32)  # Use float32 to save memory
    
    for i_start, i_end, chunk_distances in results:
        distances[i_start:i_end] = chunk_distances
    
    # Make matrix symmetric (exploit symmetry for future optimizations)
    print("Making matrix symmetric...")
    for i in range(n):
        for j in range(i + 1, n):
            # Use the average of both calculations to reduce numerical errors
            avg_dist = (distances[i, j] + distances[j, i]) / 2
            distances[i, j] = avg_dist
            distances[j, i] = avg_dist
    
    return distances


def process_symmetric_chunk_global(chunk_data):
    """Global function for multiprocessing (avoids pickle issues)"""
    pairs, feature_indices = chunk_data
    results = []
    
    # Access global features array
    global GLOBAL_FEATURES
    
    for i, j in pairs:
        distance = calculate_distance(GLOBAL_FEATURES[i], GLOBAL_FEATURES[j])
        results.append((i, j, round(distance, 5)))
    
    return results


def init_worker(features_array):
    """Initialize worker process with shared data"""
    global GLOBAL_FEATURES
    GLOBAL_FEATURES = features_array


def calculate_distances_symmetric_optimized(features, num_processes=None):
    """Optimized version using global variables to avoid pickle issues"""
    n = len(features)
    
    if num_processes is None:
        num_processes = min(os.cpu_count(), 8)
    
    print(f"Using symmetric optimization with {num_processes} processes")
    
    distances = np.zeros((n, n), dtype=np.float32)  # Use float32 to save memory
    
    # Calculate total number of unique pairs
    total_pairs = n * (n - 1) // 2
    print(f"Total unique pairs to calculate: {total_pairs}")
    
    # Create work items (i, j pairs where i < j)
    work_items = []
    for i in range(n):
        for j in range(i + 1, n):
            work_items.append((i, j))
    
    # Split work items into chunks (smaller chunks to avoid memory issues)
    chunk_size = max(100, len(work_items) // (num_processes * 8))
    chunks = []
    for i in range(0, len(work_items), chunk_size):
        chunk = work_items[i:i + chunk_size]
        # Only pass indices, not the actual features
        chunks.append((chunk, list(range(n))))
    
    print(f"Created {len(chunks)} work chunks with chunk size {chunk_size}")
    
    # Process chunks
    start_time = time.time()
    
    if num_processes == 1:
        # Sequential processing
        global GLOBAL_FEATURES
        GLOBAL_FEATURES = features
        all_results = [process_symmetric_chunk_global(chunk) for chunk in chunks]
    else:
        # Parallel processing with global variable initialization
        with Pool(processes=num_processes, initializer=init_worker, initargs=(features,)) as pool:
            all_results = pool.map(process_symmetric_chunk_global, chunks)
    
    processing_time = time.time() - start_time
    print(f"Distance calculation completed in {processing_time:.2f} seconds")
    
    # Assemble symmetric matrix
    print("Assembling symmetric matrix...")
    calculated_pairs = 0
    
    for chunk_results in all_results:
        for i, j, distance in chunk_results:
            distances[i, j] = distance
            distances[j, i] = distance
            calculated_pairs += 1
            
            # Show progress every 5% or at completion
            if calculated_pairs % max(1, total_pairs // 20) == 0 or calculated_pairs == total_pairs:
                progress = (calculated_pairs / total_pairs) * 100
                print(f"Distance calculation: {progress:.1f}% ({calculated_pairs:,}/{total_pairs:,} pairs)", end='\r')
    
    # Set diagonal to 0
    np.fill_diagonal(distances, 0)
    
    print()  # New line after progress completion
    return distances


PROGRESS_COUNTER = 0
PROGRESS_LOCK = None

def process_row_chunk_global(chunk_info):
    """Global function to process row chunks for multiprocessing"""
    start_row, end_row, n = chunk_info
    chunk_distances = []
    
    global GLOBAL_FEATURES, PROGRESS_COUNTER, PROGRESS_LOCK
    
    for i in range(start_row, end_row):
        row_distances = np.zeros(n, dtype=np.float32)  # Use float32 to save memory
        F1 = GLOBAL_FEATURES[i]
        for j in range(i + 1, n):  # Only calculate upper triangle
            distance = calculate_distance(F1, GLOBAL_FEATURES[j])
            row_distances[j] = round(distance, 5)
        chunk_distances.append((i, row_distances))
        
        # Update progress counter
        PROGRESS_COUNTER += 1
        if PROGRESS_COUNTER % 50 == 0 or PROGRESS_COUNTER == n:
            progress = (PROGRESS_COUNTER / n) * 100
            print(f"Progress: {progress:.1f}% ({PROGRESS_COUNTER}/{n} rows)", end='\r')
    
    return chunk_distances


def calculate_distances_chunked_simple(features, num_processes=None):
    """Simpler chunked approach that processes row by row"""
    n = len(features)
    
    if num_processes is None:
        num_processes = min(os.cpu_count(), 4)  # Reduced for stability
    
    # Reset progress counter
    global PROGRESS_COUNTER
    PROGRESS_COUNTER = 0
    
    distances = np.zeros((n, n), dtype=np.float32)  # Use float32 to save memory
    
    # Process in row chunks
    chunk_size = max(1, n // (num_processes * 2))
    row_chunks = []
    
    for i in range(0, n, chunk_size):
        end = min(i + chunk_size, n)
        row_chunks.append((i, end, n))  # Include n for the global function
    
    start_time = time.time()
    
    if num_processes == 1:
        global GLOBAL_FEATURES
        GLOBAL_FEATURES = features
        all_results = [process_row_chunk_global(chunk) for chunk in row_chunks]
    else:
        with Pool(processes=num_processes, initializer=init_worker, initargs=(features,)) as pool:
            all_results = pool.map(process_row_chunk_global, row_chunks)
    
    processing_time = time.time() - start_time
    print()  # New line after progress
    
    # Assemble matrix and make it symmetric
    for chunk_results in all_results:
        for i, row_distances in chunk_results:
            distances[i] = row_distances
            # Mirror to lower triangle
            for j in range(i + 1, n):
                distances[j, i] = distances[i, j]
    return distances


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Compute pairwise distances of protein sequences from phase correlation (Optimized)')
    parser.add_argument('--input', '-i', required=True, help='Input FASTA file')
    parser.add_argument('--output', '-o', required=True, help='Output CSV file path')
    parser.add_argument('--dim', '-d', type=int, default=1024, help='FFT dimension (default: 1024)')
    parser.add_argument('--processes', '-p', type=int, default=None, help='Number of processes (default: auto)')
    parser.add_argument('--method', '-m', choices=['symmetric'], default='symmetric', 
                       help='Optimization method: symmetric (always optimal)')
    args = parser.parse_args()

    dataset_path = args.input
    output_csv = args.output
    fft_dim = args.dim
    num_processes = args.processes
    method = args.method
    
    print(f"Using symmetric matrix optimization (50% computation reduction)")
    print(f"FFT dimension: {fft_dim}")
    
    embedder = ProteinEmbedder(n=fft_dim)

    ids = []
    features = []

    print("Encoding sequences...")
    start_time = time.time()
    
    for i, record in enumerate(SeqIO.parse(dataset_path, 'fasta'), 1):
        v = embedder.encode(str(record.seq))
        ids.append(str(record.id))
        features.append(v)
        
        if i % 100 == 0 or i == 1:
            print(f"Encoded {i} sequences...", end='\r')
    
    encoding_time = time.time() - start_time
    print(f"\nCompleted encoding {len(ids)} sequences in {encoding_time:.2f} seconds")

    n = len(ids)
    print(f'Number of sequences: {n}')
    print(f'Total pairwise calculations needed: {n * n:,}')
    print(f'Unique pairs (with symmetry): {n * (n - 1) // 2:,}')

    print("\nCalculating pairwise distances...")
    start_time = time.time()
    
    # Always use symmetric optimization for best performance
    distances = calculate_distances_symmetric_optimized(features, num_processes)
    
    calculation_time = time.time() - start_time
    print(f"\nCompleted calculating all pairwise distances in {calculation_time:.2f} seconds")

    print("Saving results to CSV...")
    save_start = time.time()
    
    df = pd.DataFrame(distances, index=ids, columns=ids)
    df.to_csv(output_csv, index_label='ID')
    
    save_time = time.time() - save_start
    print(f"Results saved to: {output_csv} in {save_time:.2f} seconds")
    
    total_time = encoding_time + calculation_time + save_time
    print(f"\nTotal runtime: {total_time:.2f} seconds")
    print(f"  - Encoding: {encoding_time:.2f}s ({encoding_time/total_time*100:.1f}%)")
    print(f"  - Calculation: {calculation_time:.2f}s ({calculation_time/total_time*100:.1f}%)")
    print(f"  - Saving: {save_time:.2f}s ({save_time/total_time*100:.1f}%)")

# Usage examples:
# python3 get_phase_dist_mat_optimized.py -d 512 -i phosphatase.fa -o score_matrix_optimized.csv -m symmetric
# python3 get_phase_dist_mat_optimized.py -d 512 -i phosphatase.fa -o score_matrix_optimized.csv -m vectorized -p 4