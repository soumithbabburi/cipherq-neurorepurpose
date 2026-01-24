#!/usr/bin/env python3
"""
PDB STRUCTURE HANDLER - DATABASE INTEGRATED VERSION
Downloads and processes protein structures and drug SMILES 
by querying the PostgreSQL database instead of static JSON files.
"""

import os
import requests
import io
import logging
import psycopg2
from typing import Dict, List, Any, Optional, Tuple
from Bio.PDB import PDBList

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class PDBStructureHandler:
    """
    Handle PDB structure downloading and processing for molecular docking.
    Queries the PostgreSQL database for target PDB IDs and drug SMILES.
    """
    
    def __init__(self):
        self.pdb_base_url = "https://files.rcsb.org/download"
        self.pdb_info_url = "https://data.rcsb.org/rest/v1/core/entry"
        self.structure_dir = "data/structures"
        
        # Ensure download directory exists
        if not os.path.exists(self.structure_dir):
            os.makedirs(self.structure_dir)

    def _get_db_connection(self):
        """Standard connection to the cipherq_repurpose database."""
        try:
            return psycopg2.connect(
                host="localhost",
                database="cipherq_repurpose",
                user="babburisoumith",
                password=""
            )
        except Exception as e:
            logger.error(f"❌ PDB Handler: Database connection failed: {e}")
            return None

    def get_target_structure(self, gene_symbol: str) -> Optional[str]:
        """
        1. Queries the DATABASE for the PDB ID associated with the gene.
        2. Downloads the PDB structure from RCSB if not found locally.
        """
        conn = self._get_db_connection()
        if not conn:
            return None
            
        pdb_id = None
        try:
            cur = conn.cursor()
            # Fetch PDB ID and Name from proteins table
            cur.execute("""
                SELECT pdb_id, name 
                FROM proteins 
                WHERE gene_symbol = %s 
                LIMIT 1
            """, (gene_symbol.upper(),))
            
            row = cur.fetchone()
            if row:
                pdb_id = row[0]
                logger.info(f"✅ Found PDB ID {pdb_id} for gene {gene_symbol} in DB.")
            cur.close()
        except Exception as e:
            logger.error(f"❌ Error fetching PDB ID from DB: {e}")
        finally:
            conn.close()

        if not pdb_id:
            logger.warning(f"⚠️ No PDB ID mapping found for gene: {gene_symbol}")
            return None

        # Download the actual structure file
        return self.download_pdb(pdb_id)

    def download_pdb(self, pdb_id: str) -> Optional[str]:
        """Downloads a PDB file from RCSB and returns the local file path."""
        pdb_id = pdb_id.upper()
        local_path = os.path.join(self.structure_dir, f"{pdb_id}.pdb")
        
        if os.path.exists(local_path):
            logger.info(f"📂 Using cached PDB file: {local_path}")
            return local_path
            
        try:
            logger.info(f"🌐 Downloading PDB {pdb_id} from RCSB...")
            pdbl = PDBList()
            # This downloads to the directory and returns the path (usually .ent format)
            file_path = pdbl.retrieve_pdb_file(pdb_id, pdir=self.structure_dir, file_format='pdb')
            
            # PDBList often names files 'pdbXXXX.ent', we rename to 'XXXX.pdb' for consistency
            if os.path.exists(file_path):
                os.rename(file_path, local_path)
                return local_path
        except Exception as e:
            logger.error(f"❌ Failed to download PDB {pdb_id}: {e}")
            
        return None

    def get_drug_smiles(self, drug_name: str) -> Optional[str]:
        """Queries the database to get the SMILES string for a specific drug."""
        conn = self._get_db_connection()
        if not conn:
            return None
            
        smiles = None
        try:
            cur = conn.cursor()
            # Use ILIKE for case-insensitive matching
            cur.execute("SELECT smiles FROM drugs WHERE name ILIKE %s LIMIT 1", (drug_name,))
            row = cur.fetchone()
            if row:
                smiles = row[0]
            cur.close()
        except Exception as e:
            logger.error(f"❌ Error fetching SMILES for {drug_name}: {e}")
        finally:
            conn.close()
            
        return smiles

    def prepare_docking_pair(self, drug_name: str, gene_symbol: str) -> Dict[str, Any]:
        """
        Coordinates the full preparation: 
        Fetches structure file path + Drug SMILES.
        """
        structure_path = self.get_target_structure(gene_symbol)
        smiles = self.get_drug_smiles(drug_name)
        
        if not structure_path or not smiles:
            return {"status": "error", "message": "Missing structure or SMILES data"}
            
        return {
            "status": "success",
            "drug": drug_name,
            "gene": gene_symbol,
            "smiles": smiles,
            "pdb_path": structure_path
        }

def get_pdb_handler() -> PDBStructureHandler:
    """Utility function to get an instance of the handler."""
    return PDBStructureHandler()

# Example usage for testing
if __name__ == "__main__":
    handler = get_pdb_handler()
    
    # Test case: Query Alzheimer's target (ensure ACHE or similar is in your DB)
    result = handler.prepare_docking_pair("Donepezil", "ACHE")
    print(f"Preparation result: {result}")