"""
INTELLIGENT DISEASE CONNECTION FILTER - PATHWAY-BASED APPROACH
This replaces the overly-strict direct gene matching approach

Key improvements:
1. Multi-level connection scoring (direct, pathway, biological process)
2. Uses known disease-relevant pathways
3. Considers drug mechanism of action
4. Provides confidence scores instead of binary filtering
"""

import logging

logger = logging.getLogger(__name__)

# ============================================================================
# DISEASE GENE SETS - Curated from literature
# ============================================================================

DISEASE_GENES = {
    'Alzheimers': {
        'core': {'APP', 'APOE', 'PSEN1', 'PSEN2', 'MAPT', 'BACE1'},
        'associated': {'TREM2', 'CD33', 'CLU', 'CR1', 'ABCA7', 'BIN1', 'PICALM', 
                      'MS4A6A', 'EPHA1', 'CD2AP', 'SORL1', 'PTK2B', 'SLC24A4',
                      'CELF1', 'FERMT2', 'HLA-DRB5', 'INPP5D', 'MEF2C', 'NME8',
                      'ZCWPW1', 'CASS4'},
        'pathways': {
            # Metabolic pathways (CRITICAL for diabetic drug connections)
            'metabolism': {'PRKAA1', 'PRKAA2', 'PRKAB1', 'PRKAB2', 'PRKAG1', 'PRKAG2', 'PRKAG3',  # AMPK pathway
                          'INS', 'INSR', 'IRS1', 'IRS2', 'PIK3CA', 'AKT1', 'AKT2', 'MTOR',  # Insulin signaling
                          'GSK3A', 'GSK3B',  # Glycogen synthase kinase (tau phosphorylation)
                          'PPARG', 'PPARA', 'PPARD'},  # Nuclear receptors
            
            # Mitochondrial function (Metformin's primary targets)
            'mitochondrial': {'NDUFA', 'NDUFB', 'NDUFC', 'NDUFV', 'NDUFS',  # Complex I (prefix matching)
                            'MT-ND1', 'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND5', 'MT-ND6',
                            'COX', 'CYTB', 'ATP5', 'UQCR'},  # Other complexes
            
            # Inflammation pathways
            'inflammation': {'TNF', 'IL1B', 'IL6', 'IL10', 'NFKB1', 'RELA', 'PTGS2',
                           'TREM2', 'CD33', 'TYROBP', 'CX3CR1', 'CCL2'},
            
            # Autophagy & protein clearance
            'autophagy': {'BECN1', 'ATG5', 'ATG7', 'MAP1LC3B', 'SQSTM1', 'MTOR',
                         'ULK1', 'LAMP1', 'LAMP2', 'TFEB'},
            
            # Synaptic function
            'synaptic': {'SNAP25', 'SYT1', 'SYN1', 'DLG4', 'GRIN1', 'GRIN2A', 'GRIN2B',
                        'GRIA1', 'GRIA2', 'SLC1A2', 'SLC1A3'},
            
            # Oxidative stress
            'oxidative_stress': {'SOD1', 'SOD2', 'CAT', 'GPX1', 'HMOX1', 'NQO1', 'NFE2L2'},
            
            # Lipid metabolism (APOE pathway)
            'lipid': {'APOE', 'APOA1', 'APOB', 'LDLR', 'ABCA1', 'ABCG1', 'LCAT', 'CETP',
                     'SCARB1', 'LPL', 'LIPC'}
        }
    },
    
    'Parkinsons': {
        'core': {'SNCA', 'LRRK2', 'PARK7', 'PINK1', 'PRKN', 'GBA'},
        'associated': {'MAPT', 'APOE', 'GCH1', 'VPS35', 'DNAJC6', 'SYNJ1', 'CHCHD2'},
        'pathways': {
            'dopamine': {'TH', 'DDC', 'DBH', 'DRD1', 'DRD2', 'DRD3', 'DRD4', 'DRD5',
                        'SLC6A3', 'COMT', 'MAO', 'MAOA', 'MAOB'},
            'mitochondrial': {'PINK1', 'PRKN', 'NDUFA', 'NDUFB', 'COX', 'ATP5'},
            'lysosomal': {'GBA', 'LAMP2', 'ATP13A2', 'LRRK2', 'VPS35'},
            'inflammation': {'TNF', 'IL1B', 'IL6', 'NFKB1'},
            'metabolism': {'PRKAA1', 'PRKAB1', 'PPARG', 'PPARA'}
        }
    },
    
    'Diabetic': {
        'core': {'INS', 'INSR', 'IRS1', 'IRS2', 'PPARG', 'GCK', 'KCNJ11', 'ABCC8'},
        'associated': {'TCF7L2', 'SLC30A8', 'HHEX', 'CDKAL1', 'IGF2BP2', 'CDKN2A', 'CDKN2B'},
        'pathways': {
            'insulin_signaling': {'INS', 'INSR', 'IRS1', 'IRS2', 'PIK3CA', 'PIK3R1', 'AKT1', 'AKT2',
                                 'PTEN', 'FOXO1', 'GSK3B', 'MAPK1', 'MAPK3'},
            'glucose_metabolism': {'GCK', 'PYGL', 'PFKM', 'ALDOA', 'ENO1', 'PKM', 'G6PC', 'PCK1'},
            'incretin': {'GLP1R', 'GIPR', 'DPP4', 'GCG', 'GIP'},
            'beta_cell': {'KCNJ11', 'ABCC8', 'CACNA1D', 'CACNA1C', 'SLC2A2'},
            'metabolism': {'PRKAA1', 'PRKAB1', 'PPARG', 'PPARA', 'SREBF1', 'MLXIPL'}
        }
    }
}

# ============================================================================
# PATHWAY RELEVANCE SCORING
# ============================================================================

def calculate_pathway_score(drug_targets, disease_name):
    """
    Calculate how relevant a drug's targets are to disease pathways
    Returns: (score, matched_pathways, target_details)
    """
    if disease_name not in DISEASE_GENES:
        return 0.0, [], []
    
    disease_data = DISEASE_GENES[disease_name]
    score = 0.0
    matched_pathways = []
    target_details = []
    
    # 1. Check direct gene matches (highest weight)
    core_matches = [t for t in drug_targets if t in disease_data['core']]
    if core_matches:
        score += len(core_matches) * 10.0
        target_details.append(f"Core genes: {core_matches}")
    
    associated_matches = [t for t in drug_targets if t in disease_data['associated']]
    if associated_matches:
        score += len(associated_matches) * 5.0
        target_details.append(f"Associated genes: {associated_matches}")
    
    # 2. Check pathway matches (medium weight)
    for pathway_name, pathway_genes in disease_data['pathways'].items():
        pathway_hits = []
        
        for target in drug_targets:
            # Exact match
            if target in pathway_genes:
                pathway_hits.append(target)
            # Prefix match (for gene families like NDUFA1, NDUFA2, etc.)
            elif any(target.startswith(prefix) for prefix in pathway_genes if len(prefix) > 3):
                pathway_hits.append(target)
        
        if pathway_hits:
            pathway_score = len(pathway_hits) * 2.0
            score += pathway_score
            matched_pathways.append(pathway_name)
            target_details.append(f"{pathway_name} pathway: {pathway_hits[:5]}")
    
    return score, matched_pathways, target_details

# ============================================================================
# CROSS-DISEASE CONNECTION SCORING
# ============================================================================

def calculate_cross_disease_score(source_category, target_disease, drug_targets):
    """
    Calculate repurposing potential based on known biological connections
    Example: Diabetic drugs → Alzheimer's (through metabolism, inflammation)
    """
    
    # Define known cross-disease connections
    CROSS_DISEASE_CONNECTIONS = {
        ('Diabetic', 'Alzheimers'): {
            'pathways': ['metabolism', 'mitochondrial', 'inflammation', 'autophagy', 'oxidative_stress'],
            'rationale': 'Insulin signaling, glucose metabolism, and mitochondrial function link diabetes and neurodegeneration',
            'weight': 1.5  # Boost for well-established connection
        },
        ('Cardiovascular', 'Alzheimers'): {
            'pathways': ['inflammation', 'lipid', 'oxidative_stress'],
            'rationale': 'Vascular health affects brain perfusion and Aβ clearance',
            'weight': 1.3
        },
        ('Parkinsons', 'Alzheimers'): {
            'pathways': ['mitochondrial', 'autophagy', 'inflammation', 'oxidative_stress'],
            'rationale': 'Shared protein aggregation and mitochondrial dysfunction mechanisms',
            'weight': 1.4
        }
    }
    
    connection_key = (source_category, target_disease)
    if connection_key not in CROSS_DISEASE_CONNECTIONS:
        return 0.0, []
    
    connection = CROSS_DISEASE_CONNECTIONS[connection_key]
    
    # Calculate base pathway score
    base_score, matched_pathways, _ = calculate_pathway_score(drug_targets, target_disease)
    
    # Apply weight if relevant pathways are matched
    relevant_pathways = [p for p in matched_pathways if p in connection['pathways']]
    if relevant_pathways:
        boosted_score = base_score * connection['weight']
        return boosted_score, relevant_pathways
    
    return 0.0, []

# ============================================================================
# MAIN FILTERING FUNCTION
# ============================================================================

def filter_drugs_by_disease_connection(drugs, target_disease, source_category=None, 
                                       min_score=1.0, return_all_with_scores=False):
    """
    Filter drugs based on biological connection to target disease
    
    Args:
        drugs: List of drug dicts with 'drug_name' and 'targets' (list of gene symbols)
        target_disease: Disease name (e.g., 'Alzheimers', 'Parkinsons')
        source_category: Original category if doing cross-disease repurposing (e.g., 'Diabetic')
        min_score: Minimum connection score to pass filter (default 1.0)
        return_all_with_scores: If True, return all drugs with their scores instead of filtering
    
    Returns:
        List of drugs (with added 'connection_score' field)
    """
    
    results = []
    
    for drug in drugs:
        drug_name = drug.get('drug_name', 'Unknown')
        drug_targets = drug.get('targets', [])
        
        # Extract gene symbols if targets are dicts
        if drug_targets and isinstance(drug_targets[0], dict):
            drug_targets = [t.get('gene_symbol', t.get('target', '')) for t in drug_targets]
        
        # Calculate pathway-based score
        pathway_score, matched_pathways, target_details = calculate_pathway_score(
            drug_targets, target_disease
        )
        
        # Calculate cross-disease bonus if doing repurposing
        cross_disease_score = 0.0
        cross_disease_pathways = []
        if source_category and source_category != target_disease:
            cross_disease_score, cross_disease_pathways = calculate_cross_disease_score(
                source_category, target_disease, drug_targets
            )
        
        # Total score
        total_score = pathway_score + cross_disease_score
        
        # Add metadata to drug
        drug_with_score = drug.copy()
        drug_with_score['connection_score'] = total_score
        drug_with_score['pathway_score'] = pathway_score
        drug_with_score['cross_disease_score'] = cross_disease_score
        drug_with_score['matched_pathways'] = matched_pathways
        drug_with_score['cross_disease_pathways'] = cross_disease_pathways
        drug_with_score['target_details'] = target_details
        
        # Filter or keep all
        if return_all_with_scores or total_score >= min_score:
            results.append(drug_with_score)
    
    # Sort by score
    results.sort(key=lambda x: x['connection_score'], reverse=True)
    
    # Log results
    if results:
        logger.info(f"✅ FILTER RESULTS: {len(results)}/{len(drugs)} drugs passed disease connection filter")
        if len(results) > 0:
            top_drug = results[0]
            logger.info(f"   Top drug: {top_drug.get('drug_name')} (score: {top_drug['connection_score']:.1f})")
            logger.info(f"   Pathways: {top_drug['matched_pathways']}")
    else:
        logger.warning(f"⚠️  FILTER RESULTS: 0/{len(drugs)} drugs passed filter (min_score={min_score})")
        logger.warning(f"   Consider lowering min_score or checking target data")
    
    return results

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def get_connection_explanation(drug, target_disease):
    """
    Generate human-readable explanation of why a drug is connected to a disease
    """
    explanation = []
    
    if drug.get('pathway_score', 0) > 0:
        pathways = drug.get('matched_pathways', [])
        explanation.append(f"Targets {len(pathways)} relevant pathways: {', '.join(pathways)}")
    
    if drug.get('cross_disease_score', 0) > 0:
        cross_pathways = drug.get('cross_disease_pathways', [])
        explanation.append(f"Cross-disease connection through: {', '.join(cross_pathways)}")
    
    target_details = drug.get('target_details', [])
    if target_details:
        explanation.extend(target_details)
    
    return "\n".join(explanation)

def analyze_category_for_disease(drug_list, target_disease, source_category):
    """
    Analyze an entire category of drugs for repurposing potential
    Returns summary statistics
    """
    scored_drugs = filter_drugs_by_disease_connection(
        drug_list, 
        target_disease, 
        source_category,
        return_all_with_scores=True
    )
    
    high_potential = [d for d in scored_drugs if d['connection_score'] >= 5.0]
    medium_potential = [d for d in scored_drugs if 1.0 <= d['connection_score'] < 5.0]
    low_potential = [d for d in scored_drugs if 0.1 <= d['connection_score'] < 1.0]
    
    return {
        'total_drugs': len(drug_list),
        'high_potential': len(high_potential),
        'medium_potential': len(medium_potential),
        'low_potential': len(low_potential),
        'no_connection': len(scored_drugs) - len(high_potential) - len(medium_potential) - len(low_potential),
        'top_drugs': [{'name': d['drug_name'], 'score': d['connection_score']} 
                     for d in high_potential[:10]]
    }

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    # Example: Analyze diabetic drugs for Alzheimer's repurposing
    
    example_drugs = [
        {
            'drug_name': 'Metformin',
            'targets': ['PRKAB1', 'PRKAA1', 'NDUFA11', 'NDUFB9', 'NDUFC1', 'MT-ND4']
        },
        {
            'drug_name': 'Pioglitazone',
            'targets': ['PPARG', 'PPARA', 'CYP2C8', 'SLCO1B3']
        },
        {
            'drug_name': 'Glyburide',
            'targets': ['ABCC8', 'KCNJ11', 'ABCB1']
        }
    ]
    
    print("=== ANALYZING DIABETIC DRUGS FOR ALZHEIMER'S REPURPOSING ===\n")
    
    filtered = filter_drugs_by_disease_connection(
        example_drugs,
        target_disease='Alzheimers',
        source_category='Diabetic',
        min_score=1.0
    )
    
    for drug in filtered:
        print(f"\n✅ {drug['drug_name']}")
        print(f"   Total Score: {drug['connection_score']:.1f}")
        print(f"   Pathway Score: {drug['pathway_score']:.1f}")
        print(f"   Cross-disease Score: {drug['cross_disease_score']:.1f}")
        print(f"   Matched Pathways: {drug['matched_pathways']}")
        print(f"\n   Explanation:")
        print(f"   {get_connection_explanation(drug, 'Alzheimers')}")
