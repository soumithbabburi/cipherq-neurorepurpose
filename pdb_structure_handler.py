"""
PDB Structure Handler - Uses Database Utils
No more "connection already closed" errors!
"""
import logging
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

logger = logging.getLogger(__name__)

class PDBStructureHandler:
    """Handle PDB structures and fetch SMILES from database"""
    
    def __init__(self):
        logger.info("PDB Structure Handler initialized (using database_utils)")
    
    def get_smiles_for_drug(self, drug_name: str) -> str:
        """Get SMILES from database using connection pool"""
        try:
            from database_utils import execute_query
            
            results = execute_query(
                "SELECT smiles FROM drugs WHERE name = %s LIMIT 1",
                (drug_name,)
            )
            
            if results and results[0].get('smiles'):
                smiles = results[0]['smiles']
                logger.info(f"✅ Found SMILES for {drug_name}")
                return smiles
            else:
                logger.warning(f"No SMILES found for {drug_name}")
                return None
                
        except Exception as e:
            logger.error(f"❌ Error fetching SMILES for {drug_name}: {e}")
            return None
    
    def get_protein_pdb(self, protein_name: str) -> str:
        """Get PDB ID for protein"""
        try:
            from database_utils import execute_query
            
            results = execute_query(
                "SELECT preferred_pdb, pdb_ids FROM proteins WHERE name = %s OR gene_symbol = %s LIMIT 1",
                (protein_name, protein_name)
            )
            
            if results:
                pdb_id = results[0].get('preferred_pdb')
                if not pdb_id and results[0].get('pdb_ids'):
                    pdb_ids = results[0].get('pdb_ids')
                    if pdb_ids and len(pdb_ids) > 0:
                        pdb_id = pdb_ids[0]
                
                return pdb_id
            
            return None
            
        except Exception as e:
            logger.error(f"Error fetching PDB for {protein_name}: {e}")
            return None
