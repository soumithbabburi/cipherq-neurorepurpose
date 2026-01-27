"""
JSON-Based Database Queries
Reads from JSON files instead of PostgreSQL
Much simpler - no connection issues!
"""
import json
import logging
from typing import List, Dict, Optional

logger = logging.getLogger(__name__)

# Load JSON data once at startup
_drugs = None
_proteins = None
_interactions = None

def load_json_data():
    """Load JSON files into memory"""
    global _drugs, _proteins, _interactions
    
    if _drugs is not None:
        return  # Already loaded
    
    try:
        print("Loading JSON data files...")
        
        with open('drugs.json', 'r') as f:
            _drugs = json.load(f)
        print(f"✅ Loaded {len(_drugs):,} drugs")
        
        with open('proteins.json', 'r') as f:
            _proteins = json.load(f)
        print(f"✅ Loaded {len(_proteins):,} proteins")
        
        with open('interactions.json', 'r') as f:
            _interactions = json.load(f)
        print(f"✅ Loaded {len(_interactions):,} drug-protein mappings")
        
    except FileNotFoundError as e:
        logger.error(f"JSON files not found: {e}")
        logger.error("Run extract_to_json.py first to create the JSON files!")
        _drugs = {}
        _proteins = {}
        _interactions = {}
    except Exception as e:
        logger.error(f"Failed to load JSON: {e}")
        _drugs = {}
        _proteins = {}
        _interactions = {}

# Load data on import
load_json_data()


def get_drug_targets(drug_name: str, limit: int = 10) -> List[Dict]:
    """Get targets for a drug from JSON"""
    
    # Case-insensitive search
    drug_name_lower = drug_name.lower()
    
    for name, targets in _interactions.items():
        if name.lower() == drug_name_lower:
            result = []
            for target in targets[:limit]:
                result.append({
                    'gene_symbol': target['gene_symbol'],
                    'protein_name': target['protein_name'],
                    'confidence_score': target['confidence_score'],
                    'binding_affinity': target['binding_affinity']
                })
            
            if result:
                logger.info(f"✅ {drug_name} targets: {[r['gene_symbol'] for r in result]}")
            
            return result
    
    return []


def get_drug_by_name(drug_name: str) -> Optional[Dict]:
    """Get drug info from JSON"""
    
    drug_name_lower = drug_name.lower()
    
    for name, drug_info in _drugs.items():
        if name.lower() == drug_name_lower:
            return drug_info
    
    return None


def get_drugs_by_category(category: str, limit: int = 50) -> List[Dict]:
    """
    Get drugs by category
    Since JSON doesn't have categories, return drugs that have interactions
    """
    
    # Return drugs that have targets (these are likely good candidates)
    drugs_with_targets = []
    
    for drug_name in _interactions.keys():
        if drug_name in _drugs:
            drug_info = _drugs[drug_name].copy()
            drug_info['name'] = drug_name
            drugs_with_targets.append(drug_info)
            
            if len(drugs_with_targets) >= limit:
                break
    
    return drugs_with_targets


def search_drugs_by_query(query: str, limit: int = 20) -> List[Dict]:
    """Search drugs by query"""
    
    query_lower = query.lower()
    results = []
    
    for drug_name in _interactions.keys():
        if query_lower in drug_name.lower():
            if drug_name in _drugs:
                drug_info = _drugs[drug_name].copy()
                drug_info['name'] = drug_name
                results.append(drug_info)
                
                if len(results) >= limit:
                    break
    
    return results


def get_interaction_count(drug_name: str) -> int:
    """Get number of protein interactions for a drug"""
    
    drug_name_lower = drug_name.lower()
    
    for name, targets in _interactions.items():
        if name.lower() == drug_name_lower:
            return len(targets)
    
    return 0


def get_disease_targets(disease_name: str, limit: int = 10) -> List[str]:
    """Get common protein targets for a disease"""
    disease_targets = {
        'Alzheimer': ['ACHE', 'BCHE', 'NMDA', 'APP', 'MAPT'],
        'Diabetes': ['PPARG', 'DPP4', 'SGLT2', 'GLP1R', 'AMPK'],
        'Cardiovascular': ['ACE', 'ADRB1', 'HMGCR', 'AGTR1'],
        'Cancer': ['ABL1', 'EGFR', 'ERBB2', 'KIT'],
    }
    
    for key, targets in disease_targets.items():
        if key.lower() in disease_name.lower():
            return targets[:limit]
    
    return []


# For compatibility with database_utils
def execute_query(sql, params=None, fetch=True):
    """
    Compatibility function - not used with JSON backend
    Returns empty list
    """
    logger.warning("execute_query called but using JSON backend")
    return []


__all__ = ['get_drug_targets', 'get_drug_by_name', 'get_drugs_by_category', 
           'search_drugs_by_query', 'get_interaction_count', 'get_disease_targets']


# ===== PATHWAY FUNCTIONS =====

def get_protein_pathways(gene_symbol: str) -> List[str]:
    """Get pathways for a protein"""
    if not _protein_pathways:
        return []
    
    return _protein_pathways.get(gene_symbol, [])


def get_pathway_info(pathway_id: str) -> Optional[Dict]:
    """Get pathway information"""
    if not _pathways:
        return None
    
    return _pathways.get(pathway_id)


def get_pathways_for_disease(disease_name: str) -> List[Dict]:
    """Get relevant pathways for a disease"""
    if not _pathways:
        return []
    
    disease_lower = disease_name.lower()
    relevant = []
    
    for pathway_id, pathway_info in _pathways.items():
        diseases = pathway_info.get('diseases', [])
        for disease in diseases:
            if disease.lower() in disease_lower or disease_lower in disease.lower():
                relevant.append({
                    'pathway_id': pathway_id,
                    'name': pathway_info['name'],
                    'category': pathway_info.get('category', 'Unknown'),
                    'proteins': pathway_info.get('proteins', []),
                    'relevance_score': pathway_info.get('relevance_score', 0.7)
                })
                break
    
    return sorted(relevant, key=lambda x: x['relevance_score'], reverse=True)


def get_drug_pathways(drug_name: str) -> List[Dict]:
    """Get pathways for a drug (via its protein targets)"""
    targets = get_drug_targets(drug_name)
    
    if not targets or not _pathways:
        return []
    
    pathway_set = set()
    for target in targets:
        gene_symbol = target['gene_symbol']
        pathways = get_protein_pathways(gene_symbol)
        pathway_set.update(pathways)
    
    # Get full pathway info
    pathway_list = []
    for pathway_id in pathway_set:
        pathway_info = get_pathway_info(pathway_id)
        if pathway_info:
            pathway_list.append({
                'pathway_id': pathway_id,
                'name': pathway_info['name'],
                'category': pathway_info.get('category', 'Unknown'),
                'relevance_score': pathway_info.get('relevance_score', 0.7)
            })
    
    return sorted(pathway_list, key=lambda x: x['relevance_score'], reverse=True)


__all__ = ['get_drug_targets', 'get_drug_by_name', 'get_drugs_by_category', 
           'search_drugs_by_query', 'get_interaction_count', 'get_disease_targets',
           'get_protein_pathways', 'get_pathway_info', 'get_pathways_for_disease',
           'get_drug_pathways']
