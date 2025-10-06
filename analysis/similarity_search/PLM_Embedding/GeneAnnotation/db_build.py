#!/usr/bin/env python3
"""
Main entry point for the protein vector database builder application.
"""

import logging
import sys
from embedders import get_embedder
from database import VectorDatabase
from processor import ProteinProcessor
from config import parse_args, ConfigurationError

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('protein_embedder.log')
    ]
)
logger = logging.getLogger(__name__)


def main():
    """Main function to run the protein database builder."""
    try:
        # Parse command line arguments
        args = parse_args()

        # Initialize embedder based on model name
        logger.info(f"Loading embedding model: {args.model_name}")
        embedder = get_embedder(args.model_name)
        logger.info(f'Embedding model loaded successfully: {embedder.__class__.__name__}')

        # Initialize vector database
        db = VectorDatabase(
            collection_name=args.collection,
            vector_size=embedder.vector_size
        )

        # Initialize processor with optimization parameters
        processor = ProteinProcessor(
            embedder=embedder,
            db=db,
            batch_size=args.batch_size,
            num_workers=args.num_workers,
            embedding_batch_size=args.embedding_batch_size
        )
        logger.info('Processor initialized successfully')

        # Process the FASTA file using optimized method
        stats = processor.process_fasta_file_optimized(args.fasta_path)

        logger.info("Program completed successfully")
        logger.info(f"Total sequences: {stats.total_sequences}, "
                    f"Processed: {stats.processed_sequences}, "
                    f"Failed: {stats.failed_sequences}, "
                    f"Empty: {stats.empty_sequences}")

        return 0

    except ConfigurationError as e:
        logger.error(f"Configuration error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Program failed: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
