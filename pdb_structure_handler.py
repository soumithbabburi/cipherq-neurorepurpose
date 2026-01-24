#!/usr/bin/env python3
"""
PDB STRUCTURE HANDLER - NEON DATABASE VERSION
Downloads and processes protein structures and drug SMILES 
by querying the Neon PostgreSQL database instead of local JSON files.
"""

import os
import logging
from typing import Dict, Any, Optional
from Bio.PDB import PDBList
import psycopg2
import streamlit as st

# Optional: import your Config if available
try:
    from config import Config
    CONFIG_AVAILABLE = True
except ImportError:
    CONFIG_AVAILABLE = False

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Neon connection string
NEON_CONN_STRING = (
    "postgresql://neondb_owner:npg_X3pDYGci7asM@ep-orange-dust-afdedmkv-pooler.c-2.us-west-2.aws.neon.tech/"
    "neondb?sslmode=require&channel_binding=require"
)

@st.cache_resource
def get_db_connection():
    """Create and cache a persistent connection to Neon database"""
    try:
        if CONFIG_AVAILABLE:
            db_params = Config.get_db_params()
            conn = psycopg2.connect(**db_params)
        else:
            # Connect to Neon using the full connection string
            conn = psycopg2.connect(NEON_CONN_STRING)

        logger.info(f"✅ Connected to Neon database")
        return conn
    except Exception as e:
        logger.error(f"❌ Database connection failed: {e}")
        return None

class PDBStructureHandler:
    """
    Handle PDB structure downloading and processing for molecular docking.
    Queries the Neon PostgreSQL database for target PDB IDs and drug SMILES.
    """
    
    def __init__(self):
        self.pdb_base_url = "https://files.rcsb.org/download"
        self.structure_dir = "data/structures"
        
        if not os.path.exists(self.structure_dir):
            os.makedirs(self.structure_dir)
        
        # Persistent DB connection
        self.conn = get_db_connection()
        if not self.conn:
            logger.error("❌ Could not establish Neon database connection.")

    def __del__(self):
        if self.conn:
            self.conn.close()
            logger.info("✅ Neon database connection closed.")

    def get_target_structure(self, gene_symbol: str) -> Optional[str]:
        if not self.conn:
            logger.error("❌ No active DB connection.")
            return None

        pdb_id = None
        try:
            cur = self.conn.cursor()
            cur.execute("""
                SELECT pdb_id, name 
                FROM proteins 
                WHERE gene_symbol = %s 
                LIMIT 1
            """, (gene_symbol.upper(),))
            row = cur.fetchone()
            if row:
                pdb_id = row[0]
                logger.info(f"✅ Found PDB ID {pdb_id} for gene {gene_symbol}.")
            cur.close()
        except Exception as e:
            logger.error(f"❌ Error fetching PDB ID from DB: {e}")

        if not pdb_id:
            logger.warning(f"⚠️ No PDB ID mapping found for gene: {gene_symbol}")
            return None

        return self.download_pdb(pdb_id)

    def download_pdb(self, pdb_id: str) -> Optional[str]:
        pdb_id = pdb_id.upper()
        local_path = os.path.join(self.structure_dir, f"{pdb_id}.pdb")
        
        if os.path.exists(local_path):
            logger.info(f"📂 Using cached PDB file: {local_path}")
            return local_path
            
        try:
            logger.info(f"🌐 Downloading PDB {pdb_id} from RCSB...")
            pdbl = PDBList()
            file_path = pdbl.retrieve_pdb_file(pdb_id, pdir=self.structure_dir, file_format='pdb')
            if os.path.exists(file_path):
                os.rename(file_path, local_path)
                return local_path
        except Exception as e:
            logger.error(f"❌ Failed to download PDB {pdb_id}: {e}")
            
        return None

    def get_drug_smiles(self, drug_name: str) -> Optional[str]:
        if not self.conn:
            logger.error("❌ No active DB connection.")
            return None
            
        smiles = None
        try:
            cur = self.conn.cursor()
            cur.execute("SELECT smiles FROM drugs WHERE name ILIKE %s LIMIT 1", (drug_name,))
            row = cur.fetchone()
            if row:
                smiles = row[0]
            cur.close()
        except Exception as e:
            logger.error(f"❌ Error fetching SMILES for {drug_name}: {e}")
            
        return smiles

    def prepare_docking_pair(self, drug_name: str, gene_symbol: str) -> Dict[str, Any]:
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
    return PDBStructureHandler()

# Example usage
if __name__ == "__main__":
    handler = get_pdb_handler()
    result = handler.prepare_docking_pair("Pioglitazone", "PPARG")
    print(f"Preparation result: {result}")
