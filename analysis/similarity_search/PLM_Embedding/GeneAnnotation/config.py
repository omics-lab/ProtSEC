"""
Configuration and argument parsing module for protein database builder.
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, Any
import torch
import platform


def get_device() -> torch.device:
    """Determine the appropriate device for computation."""
    if torch.cuda.is_available():
        return torch.device('cuda')
    elif torch.backends.mps.is_available() and platform.system() == 'Darwin':
        return torch.device('mps')
    return torch.device('cpu')


class ConfigurationError(Exception):
    """Raised when there's an error in the configuration."""
    pass

def parse_args() -> argparse.Namespace:
    """Parse command line arguments with input validation."""
    parser = argparse.ArgumentParser(
        description='Build protein vector database from FASTA file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--fasta_path', type=str, required=True, help='Path to input FASTA file')
    parser.add_argument('--collection', type=str, required=True, help='Name of the collection to create in the Database')
    parser.add_argument('--batch_size', type=int, default=50, help='Batch size for database uploads')
    parser.add_argument('--embedding_batch_size', type=int, default=16, help='Batch size for embedding generation')
    parser.add_argument('--num_workers', type=int, default=4, help='Number of parallel workers for processing')
    parser.add_argument('--model_name', type=str, default="esm2", choices=["prot_bert", "esm2_small", "esm2_large", "prot_t5", "openai"], help='Protein embedding model to use')

    args = parser.parse_args()

    # Validate arguments
    if not Path(args.fasta_path).exists():
        raise ConfigurationError(f"FASTA file not found: {args.fasta_path}")
    if args.batch_size < 1:
        raise ConfigurationError("Batch size must be positive")
    if args.embedding_batch_size < 1:
        raise ConfigurationError("Embedding batch size must be positive")
    if args.num_workers < 1:
        raise ConfigurationError("Number of workers must be positive")
    
    return args
