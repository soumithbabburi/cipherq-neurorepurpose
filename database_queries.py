"""
Database Queries - JSON ONLY VERSION
NO PostgreSQL - uses only JSON files
"""

import json
import logging
import os

logger = logging.getLogger(__name__)

# ===== LOAD ALL JSON FILES =====
_DRUG_INTERACTIONS = None
_PATHWAYS = None
_PROTEIN_PATHWAYS = None
_DRUGS = None
_GENES = None

def load_all_data():
    """Load all JSON files"""
    global _DRUG_INTERACTIONS, _PATHWAYS, _PROTEIN_PATHWAYS, _DRUGS, _GENES
    
    if _DRUG_INTERACTIONS is not None:
        return
    
    # Load drug interactions
    try:
        with open('drug_interactions.json', 'r') as f:
            _DRUG_INTERACTIONS = json.load(f)
        logger.info(f"✅ Loaded {len(_DRUG_INTERACTIONS)} drugs with interactions")
    except:
        _DRUG_INTERACTIONS = {}
        logger.warning("drug_interactions.json not found")
    
    # Load pathways
    try:
        with open('pathways.json', 'r') as f:
            _PATHWAYS = json.load(f)
        logger.info(f"✅ Loaded {len(_PATHWAYS)} pathways")
    except:
        _PATHWAYS = {}
        logger.warning("pathways.json not found")
    
    # Load protein-pathway mappings
    try:
        with open('protein_pathways.json', 'r') as f:
            _PROTEIN_PATHWAYS = json.load(f)
        logger.info(f"✅ Loaded {len(_PROTEIN_PATHWAYS)} protein-pathway mappings")
    except:
        _PROTEIN_PATHWAYS = {}
        logger.warning("protein_pathways.json not found")
    
    # Load drugs metadata
    try:
        with open('drugs.json', 'r') as f:
            _DRUGS = json.load(f)
        logger.info(f"✅ Loaded {len(_DRUGS)} drugs with metadata")
    except:
        _DRUGS = {}
        logger.warning("drugs.json not found")
    
    # Load official genes
    try:
        with open('genes.json', 'r') as f:
            _GENES = json.load(f)
        logger.info(f"✅ Loaded {len(_GENES)} official HGNC genes")
    except:
        _GENES = {}
        logger.warning("genes.json not found")

load_all_data()

# ===== QUERY FUNCTIONS =====

def get_drug_targets(drug_name: str):
    """Get protein targets for a drug"""
    load_all_data()
    return _DRUG_INTERACTIONS.get(drug_name.lower(), [])

def get_drug_info(drug_name: str):
    """Get drug metadata (SMILES, properties)"""
    load_all_data()
    return _DRUGS.get(drug_name.lower())

def get_gene_info(gene_symbol: str):
    """Get official gene information from HGNC"""
    load_all_data()
    return _GENES.get(gene_symbol.upper())

def search_drugs(query: str, limit: int = 50):
    """Search drugs by name"""
    load_all_data()
    query_lower = query.lower()
    
    results = []
    for drug_key, drug_data in _DRUGS.items():
        if query_lower in drug_key or query_lower in drug_data.get('name', '').lower():
            results.append(drug_data)
            if len(results) >= limit:
                break
    
    return results

def get_pathway_info(pathway_id: str):
    """Get pathway details"""
    load_all_data()
    return _PATHWAYS.get(pathway_id)

def get_gene_pathways(gene_symbol: str):
    """Get pathways for a gene"""
    load_all_data()
    return _PROTEIN_PATHWAYS.get(gene_symbol.upper(), [])

# ===== NO DATABASE - ALL JSON =====
def execute_query(*args, **kwargs):
    """Removed - use JSON files instead"""
    logger.error("execute_query called but database removed! Use JSON functions instead")
    return []

print("✅ Database queries module loaded - JSON ONLY (no PostgreSQL)")
print(f"   Loaded: drugs, genes, interactions, pathways, protein_pathways")
