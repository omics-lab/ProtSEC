"""
Module for managing the vector database connection and operations.
"""

import logging
from typing import List, Dict, Any
from pymilvus import MilvusClient, DataType
import os

logger = logging.getLogger(__name__)


class DatabaseError(Exception):
    """Raised when there's an error with database operations."""
    pass


class VectorDatabase:
    """Class to handle vector database operations."""

    def __init__(self, collection_name: str, vector_size: int, uri: str = "vector.db"):
        """Initialize the vector database connection."""
        
        self.vector_size = vector_size
        self.collection_name = collection_name

        try:
            self.client = MilvusClient(uri)
            self._initialize_collection()
        except Exception as e:
            raise DatabaseError(f"Failed to initialize database: {str(e)}")

    def _initialize_collection(self):
        """Create or recreate collection."""
        try:
            schema = MilvusClient.create_schema(
                auto_id=False,
                enable_dynamic_field=True
            )
            schema.add_field(field_name="id", datatype=DataType.INT64, is_primary=True)
            schema.add_field(field_name="vector", datatype=DataType.FLOAT_VECTOR, dim=self.vector_size)
            schema.add_field(field_name="protein_info", datatype=DataType.VARCHAR, max_length=1024)
            schema.add_field(field_name="sequence_length", datatype=DataType.INT64)

            # Prepare index parameters
            index_params = self.client.prepare_index_params()

            # Add index only for vector field (primary key is auto-indexed)
            index_params.add_index(
                field_name="vector", 
                index_type="FLAT",
                index_name="vector_index",
                metric_type="IP"
            )

            if self.client.has_collection(collection_name=self.collection_name):
                self.client.drop_collection(collection_name=self.collection_name)

            self.client.create_collection(
                collection_name=self.collection_name,
                schema=schema,
                index_params=index_params
            )

        except Exception as e:
            raise DatabaseError(f"Failed to initialize collection: {str(e)}")

    def upload_points(self, data: List[Dict]):
        """Upload points to the database."""
        if not data:
            return
        
        try:
            self.client.insert(
                collection_name=self.collection_name,
                data=data
            )
            logger.debug(f"Uploaded {len(data)} vectors to {self.collection_name}")
        except Exception as e:
            raise DatabaseError(f"Failed to upload points: {str(e)}")

    def close(self):
        """Close the database connection."""
        try:
            self.client.close()
        except Exception as e:
            logger.error(f"Error closing database connection: {str(e)}")


if __name__ == '__main__':
    client = VectorDatabase("protein_collection", 3)
    try:
        data = [
            {"id": 1, "vector": [0.1, 0.2, 0.3], "protein_info": "protein1", "sequence_length": 100},
            {"id": 2, "vector": [0.4, 0.5, 0.6], "protein_info": "protein2", "sequence_length": 200},
        ]
        client.upload_points(data)
        print("Done!")
    finally:
        client.close()
