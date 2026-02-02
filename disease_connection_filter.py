"""
ULTIMATE DISEASE CONNECTION FILTER - ALL-IN-ONE VERSION
=========================================================
Drop-in replacement for disease_connection_filter.py

Features:
✅ Pathway-based connection scoring
✅ Cross-disease repurposing support  
✅ Automatic data enrichment from interactions file
✅ Comprehensive logging and diagnostics
✅ Handles multiple data formats
✅ No external dependencies except logging and json

Usage:
    from disease_connection_filter import filter_drugs_by_disease_connection
    
    filtered = filter_drugs_by_disease_connection(
        drugs=my_drug_list,
        target_disease='Alzheimers',
        source_category='Diabetic',
        min_score=1.0
    )
"""

import logging
import json
import os

logger = logging.getLogger(__name__)

# ============================================================================
# DATA LOADING AND CACHING
# ============================================================================

_INTERACTIONS_CACHE = None

def load_interactions_database(file_path='drug_interactions.json'):
    """Load drug-target interactions from JSON file"""
    # Try multiple possible locations
    possible_paths = [
        file_path,
        'drug_interactions.json',
        'data/drug_interactions.json',
        'hetionet_drug_interactions.json',
        'data/hetionet_drug_interactions.json'
    ]
    
    for path in possible_paths:
        try:
            if os.path.exists(path):
                with open(path, 'r') as f:
                    data = json.load(f)
                    logger.info(f"✅ Loaded interactions from: {path}")
                    return data
        except Exception as e:
            continue
    
    logger.error(f"❌ Could not find interactions file in any of: {possible_paths}")
    return {}

def get_interactions_cache(force_reload=False):
    """Get or load interactions cache"""
    global _INTERACTIONS_CACHE
    if _INTERACTIONS_CACHE is None or force_reload:
        _INTERACTIONS_CACHE = load_interactions_database()
        if _INTERACTIONS_CACHE:
            logger.info(f"✅ Cached interactions for {len(_INTERACTIONS_CACHE)} drugs")
    return _INTERACTIONS_CACHE

# ============================================================================
# DISEASE PATHWAY DEFINITIONS
# ============================================================================

DISEASE_GENES = {
    'Alzheimers': {
        'core': {'APP', 'APOE', 'PSEN1', 'PSEN2', 'MAPT', 'BACE1'},
        'associated': {'TREM2', 'CD33', 'CLU', 'CR1', 'ABCA7', 'BIN1', 'PICALM', 
                      'MS4A6A', 'EPHA1', 'CD2AP', 'SORL1', 'PTK2B'},
        'pathways': {
            'metabolism': {'PRKAA1', 'PRKAA2', 'PRKAB1', 'PRKAB2', 'PRKAG1',
                          'INS', 'INSR', 'IRS1', 'IRS2', 'PIK3CA', 'AKT1', 'AKT2',
                          'GSK3A', 'GSK3B', 'PPARG', 'PPARA'},
            'mitochondrial': {'ND', 'NDUFA', 'NDUFB', 'NDUFC', 'NDUFV', 'NDUFS',
                            'MT-ND', 'COX', 'ATP5', 'UQCR'},
            'inflammation': {'TNF', 'IL1B', 'IL6', 'IL10', 'NFKB1', 'PTGS2',
                           'TREM2', 'CD33', 'CX3CR1'},
            'autophagy': {'BECN1', 'ATG5', 'ATG7', 'MAP1LC3B', 'MTOR', 'ULK1'},
            'synaptic': {'SNAP25', 'SYT1', 'SYN1', 'DLG4', 'GRIN1', 'GRIN2A'},
            'oxidative_stress': {'SOD1', 'SOD2', 'CAT', 'GPX1', 'HMOX1'},
            'lipid': {'APOE', 'APOA1', 'LDLR', 'ABCA1', 'LCAT', 'SCARB1'}
        }
    },
    'Parkinsons': {
        'core': {'SNCA', 'LRRK2', 'PARK7', 'PINK1', 'PRKN', 'GBA'},
        'associated': {'MAPT', 'APOE', 'GCH1', 'VPS35'},
        'pathways': {
            'dopamine': {'TH', 'DDC', 'DRD1', 'DRD2', 'SLC6A3', 'COMT', 'MAO'},
            'mitochondrial': {'PINK1', 'PRKN', 'NDUFA', 'NDUFB', 'COX'},
            'lysosomal': {'GBA', 'LAMP2', 'ATP13A2', 'LRRK2'},
            'inflammation': {'TNF', 'IL1B', 'IL6', 'NFKB1'},
            'metabolism': {'PRKAA1', 'PRKAB1', 'PPARG'}
        }
    },
    'Diabetic': {
        'core': {'INS', 'INSR', 'IRS1', 'IRS2', 'PPARG', 'GCK', 'KCNJ11'},
        'associated': {'TCF7L2', 'SLC30A8', 'HHEX', 'CDKAL1'},
        'pathways': {
            'insulin_signaling': {'INS', 'INSR', 'IRS1', 'IRS2', 'PIK3CA',
                                 'AKT1', 'AKT2', 'GSK3B'},
            'glucose_metabolism': {'GCK', 'PYGL', 'PFKM', 'PKM', 'G6PC'},
            'incretin': {'GLP1R', 'DPP4', 'GCG'},
            'metabolism': {'PRKAA1', 'PRKAB1', 'PPARG', 'PPARA'}
        }
    }
}

CROSS_DISEASE_CONNECTIONS = {
    ('Diabetic', 'Alzheimers'): {
        'pathways': ['metabolism', 'mitochondrial', 'inflammation', 'autophagy'],
        'weight': 1.5
    },
    ('Cardiovascular', 'Alzheimers'): {
        'pathways': ['inflammation', 'lipid', 'oxidative_stress'],
        'weight': 1.3
    },
    ('Parkinsons', 'Alzheimers'): {
        'pathways': ['mitochondrial', 'autophagy', 'inflammation'],
        'weight': 1.4
    }
}

# ============================================================================
# DATA ENRICHMENT
# ============================================================================

def extract_gene_symbols(targets):
    """Extract gene symbols from various target formats"""
    if not targets:
        return []
    
    gene_symbols = []
    for target in targets:
        if isinstance(target, str):
            gene_symbols.append(target)
        elif isinstance(target, dict):
            symbol = target.get('gene_symbol') or target.get('target') or target.get('gene')
            if symbol:
                gene_symbols.append(symbol)
    
    return gene_symbols

def normalize_drug_name(name):
    """Normalize drug name for matching"""
    if not name:
        return ""
    # Convert to lowercase, remove extra spaces and special chars
    normalized = name.lower().strip()
    normalized = normalized.replace('-', ' ')
    normalized = ' '.join(normalized.split())  # Remove extra spaces
    return normalized

def build_drug_lookup_index(interactions_db):
    """Build index for fast drug name lookup with variations"""
    lookup = {}
    
    for drug_name in interactions_db.keys():
        # Original name
        lookup[drug_name] = drug_name
        
        # Lowercase
        lookup[drug_name.lower()] = drug_name
        
        # Normalized
        normalized = normalize_drug_name(drug_name)
        lookup[normalized] = drug_name
        
        # Without spaces
        lookup[drug_name.replace(' ', '')] = drug_name
        lookup[drug_name.lower().replace(' ', '')] = drug_name
    
    return lookup

_DRUG_LOOKUP_INDEX = None

def find_drug_in_interactions(drug_name, interactions_db):
    """Find drug in interactions database with fuzzy matching"""
    global _DRUG_LOOKUP_INDEX
    
    if not drug_name:
        return None
    
    # Build lookup index if not exists
    if _DRUG_LOOKUP_INDEX is None:
        _DRUG_LOOKUP_INDEX = build_drug_lookup_index(interactions_db)
    
    # Try exact match first
    if drug_name in interactions_db:
        return drug_name
    
    # Try lookup index
    if drug_name in _DRUG_LOOKUP_INDEX:
        return _DRUG_LOOKUP_INDEX[drug_name]
    
    # Try lowercase
    if drug_name.lower() in _DRUG_LOOKUP_INDEX:
        return _DRUG_LOOKUP_INDEX[drug_name.lower()]
    
    # Try normalized
    normalized = normalize_drug_name(drug_name)
    if normalized in _DRUG_LOOKUP_INDEX:
        return _DRUG_LOOKUP_INDEX[normalized]
    
    # Try without spaces
    no_space = drug_name.replace(' ', '')
    if no_space in _DRUG_LOOKUP_INDEX:
        return _DRUG_LOOKUP_INDEX[no_space]
    
    return None

def enrich_drug_with_targets(drug_dict, interactions_db=None):
    """Add target information to drug dictionary"""
    if interactions_db is None:
        interactions_db = get_interactions_cache()
    
    # Get drug name from various possible fields
    drug_name = (drug_dict.get('drug_name') or 
                 drug_dict.get('name') or 
                 drug_dict.get('Drug') or
                 drug_dict.get('compound_name'))
    
    if not drug_name:
        logger.warning(f"⚠️ Drug missing name field: {list(drug_dict.keys())}")
        return drug_dict
    
    # Check if targets already present
    if 'targets' in drug_dict and drug_dict['targets']:
        # Ensure they're gene symbols, not dicts
        drug_dict['targets'] = extract_gene_symbols(drug_dict['targets'])
        if 'drug_name' not in drug_dict:
            drug_dict['drug_name'] = drug_name
        return drug_dict
    
    # Look up targets using fuzzy matching
    matched_name = find_drug_in_interactions(drug_name, interactions_db)
    
    if matched_name:
        drug_dict['targets'] = [t['gene_symbol'] for t in interactions_db[matched_name]]
        if matched_name != drug_name:
            logger.debug(f"Matched '{drug_name}' to '{matched_name}'")
    else:
        drug_dict['targets'] = []
        logger.debug(f"No match found for '{drug_name}'")
    
    # Ensure drug_name field
    if 'drug_name' not in drug_dict:
        drug_dict['drug_name'] = drug_name
    
    return drug_dict

# ============================================================================
# SCORING FUNCTIONS
# ============================================================================

def matches_pathway_gene(target, pathway_genes):
    """Check if target matches pathway genes (including prefix matching)"""
    if target in pathway_genes:
        return True
    
    # Check prefix matching for gene families (e.g., NDUFA1, NDUFA2 match 'NDUFA')
    for pathway_gene in pathway_genes:
        if len(pathway_gene) > 3 and target.startswith(pathway_gene):
            return True
    
    return False

def calculate_pathway_score(drug_targets, disease_name):
    """Calculate pathway-based connection score"""
    if disease_name not in DISEASE_GENES:
        return 0.0, [], []
    
    disease_data = DISEASE_GENES[disease_name]
    score = 0.0
    matched_pathways = []
    target_details = []
    
    # Direct gene matches (highest weight)
    core_matches = [t for t in drug_targets if t in disease_data['core']]
    if core_matches:
        score += len(core_matches) * 10.0
        target_details.append(f"Core: {core_matches[:3]}")
    
    associated_matches = [t for t in drug_targets if t in disease_data['associated']]
    if associated_matches:
        score += len(associated_matches) * 5.0
        target_details.append(f"Associated: {associated_matches[:3]}")
    
    # Pathway matches (medium weight)
    for pathway_name, pathway_genes in disease_data['pathways'].items():
        pathway_hits = [t for t in drug_targets if matches_pathway_gene(t, pathway_genes)]
        
        if pathway_hits:
            pathway_score = len(pathway_hits) * 2.0
            score += pathway_score
            matched_pathways.append(pathway_name)
            target_details.append(f"{pathway_name}: {pathway_hits[:3]}")
    
    return score, matched_pathways, target_details

def calculate_cross_disease_score(source_category, target_disease, drug_targets):
    """Calculate cross-disease repurposing bonus"""
    connection_key = (source_category, target_disease)
    if connection_key not in CROSS_DISEASE_CONNECTIONS:
        return 0.0, []
    
    connection = CROSS_DISEASE_CONNECTIONS[connection_key]
    base_score, matched_pathways, _ = calculate_pathway_score(drug_targets, target_disease)
    
    relevant_pathways = [p for p in matched_pathways if p in connection['pathways']]
    if relevant_pathways:
        return base_score * connection['weight'], relevant_pathways
    
    return 0.0, []

# ============================================================================
# MAIN FILTERING FUNCTION
# ============================================================================

def filter_drugs_by_disease_connection(drugs, target_disease, source_category=None, 
                                       min_score=1.0, auto_enrich=True):
    """
    Filter drugs based on biological connection to target disease
    
    Args:
        drugs: List of drug dicts (will auto-enrich with targets if needed)
        target_disease: Disease name ('Alzheimers', 'Parkinsons', etc.)
        source_category: Source category for cross-disease scoring
        min_score: Minimum connection score threshold
        auto_enrich: Automatically add targets from interactions DB
    
    Returns:
        List of drugs with connection scores
    """
    
    if not drugs:
        logger.warning("⚠️ Empty drug list provided to filter")
        return []
    
    # Auto-enrich with targets if needed
    if auto_enrich:
        interactions_db = get_interactions_cache()
        if not interactions_db:
            logger.error("❌ No interactions database loaded! Check file path.")
            return []
        
        logger.info(f"Enriching {len(drugs)} drugs with target data...")
        enriched_drugs = [enrich_drug_with_targets(d.copy(), interactions_db) for d in drugs]
    else:
        enriched_drugs = drugs
    
    results = []
    no_targets_count = 0
    drugs_with_targets = []
    
    for drug in enriched_drugs:
        drug_name = drug.get('drug_name', 'Unknown')
        drug_targets = extract_gene_symbols(drug.get('targets', []))
        
        if not drug_targets:
            no_targets_count += 1
            continue
        
        drugs_with_targets.append(drug_name)
        
        # Calculate scores
        pathway_score, matched_pathways, target_details = calculate_pathway_score(
            drug_targets, target_disease
        )
        
        cross_disease_score = 0.0
        cross_disease_pathways = []
        if source_category and source_category != target_disease:
            cross_disease_score, cross_disease_pathways = calculate_cross_disease_score(
                source_category, target_disease, drug_targets
            )
        
        total_score = pathway_score + cross_disease_score
        
        # Create result drug
        drug_with_score = drug.copy()
        drug_with_score['connection_score'] = total_score
        drug_with_score['pathway_score'] = pathway_score
        drug_with_score['cross_disease_score'] = cross_disease_score
        drug_with_score['matched_pathways'] = matched_pathways
        drug_with_score['cross_disease_pathways'] = cross_disease_pathways
        drug_with_score['target_details'] = target_details
        drug_with_score['num_targets'] = len(drug_targets)
        
        if total_score >= min_score:
            results.append(drug_with_score)
    
    # Sort by score
    results.sort(key=lambda x: x['connection_score'], reverse=True)
    
    # Detailed logging
    logger.info(f"📊 ENRICHMENT STATS:")
    logger.info(f"   Input drugs: {len(drugs)}")
    logger.info(f"   Drugs with targets: {len(drugs_with_targets)}")
    logger.info(f"   Drugs without targets: {no_targets_count}")
    
    if no_targets_count > 0 and len(drugs_with_targets) > 0:
        logger.info(f"   ✅ Sample drugs WITH targets: {drugs_with_targets[:5]}")
    
    if results:
        logger.info(f"✅ FILTER RESULTS: {len(results)}/{len(drugs)} drugs passed (min_score={min_score})")
        top = results[0]
        logger.info(f"   Top: {top.get('drug_name')} (score: {top['connection_score']:.1f}, pathways: {len(top['matched_pathways'])})")
        if len(results) > 1:
            logger.info(f"   2nd: {results[1].get('drug_name')} (score: {results[1]['connection_score']:.1f})")
    else:
        logger.warning(f"⚠️ FILTER RESULTS: 0/{len(drugs)} drugs passed filter (min_score={min_score})")
        if no_targets_count > 0:
            logger.warning(f"   {no_targets_count} drugs had no target data")
        if no_targets_count == len(drugs):
            logger.error(f"   🔴 CRITICAL: NO drugs have target data!")
            logger.error(f"   Check: (1) Is drug_interactions.json in repo root? (2) Are drug names correct?")
        logger.warning(f"   Consider: (1) lowering min_score, (2) checking interaction data, (3) using auto_enrich=True")
    
    return results

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def get_connection_explanation(drug, target_disease):
    """Generate human-readable connection explanation"""
    parts = []
    
    if drug.get('pathway_score', 0) > 0:
        pathways = drug.get('matched_pathways', [])
        parts.append(f"Targets {len(pathways)} pathways: {', '.join(pathways[:3])}")
    
    if drug.get('cross_disease_score', 0) > 0:
        cross = drug.get('cross_disease_pathways', [])
        parts.append(f"Cross-disease boost via: {', '.join(cross)}")
    
    details = drug.get('target_details', [])
    if details:
        parts.extend(details[:3])
    
    return " | ".join(parts)

def summarize_filtering_results(results, source_category=None):
    """Generate summary statistics of filtering results"""
    if not results:
        return {
            'total': 0,
            'high_score': 0,
            'medium_score': 0,
            'low_score': 0
        }
    
    high = [d for d in results if d['connection_score'] >= 5.0]
    medium = [d for d in results if 1.0 <= d['connection_score'] < 5.0]
    low = [d for d in results if d['connection_score'] < 1.0]
    
    summary = {
        'total': len(results),
        'high_score': len(high),
        'medium_score': len(medium),
        'low_score': len(low),
        'top_drugs': [
            {
                'name': d['drug_name'],
                'score': d['connection_score'],
                'pathways': d['matched_pathways']
            }
            for d in results[:10]
        ]
    }
    
    return summary

# ============================================================================
# BACKWARD COMPATIBILITY
# ============================================================================

# Keep old function signature for compatibility
def filter_drugs_for_disease(drugs, disease_name, min_score=1.0):
    """Legacy function signature"""
    return filter_drugs_by_disease_connection(
        drugs,
        target_disease=disease_name,
        min_score=min_score
    )

# ============================================================================
# TESTING
# ============================================================================

if __name__ == "__main__":
    print("="*70)
    print("TESTING DISEASE CONNECTION FILTER")
    print("="*70)
    
    # Test with sample data
    test_drugs = [
        {'drug_name': 'Metformin', 'category': 'Diabetic'},
        {'drug_name': 'Pioglitazone', 'category': 'Diabetic'},
        {'drug_name': 'Glyburide', 'category': 'Diabetic'}
    ]
    
    print("\nTest 1: Diabetic → Alzheimers repurposing")
    filtered = filter_drugs_by_disease_connection(
        test_drugs,
        target_disease='Alzheimers',
        source_category='Diabetic',
        min_score=1.0
    )
    
    for i, drug in enumerate(filtered, 1):
        print(f"\n{i}. {drug['drug_name']}")
        print(f"   Score: {drug['connection_score']:.1f}")
        print(f"   Pathways: {', '.join(drug['matched_pathways'])}")
        print(f"   Explanation: {get_connection_explanation(drug, 'Alzheimers')}")
    
    print("\n" + "="*70)
    summary = summarize_filtering_results(filtered)
    print(f"SUMMARY: {summary['total']} drugs passed filter")
    print(f"  High (≥5.0): {summary['high_score']}")
    print(f"  Medium (1-5): {summary['medium_score']}")
    print("="*70)
