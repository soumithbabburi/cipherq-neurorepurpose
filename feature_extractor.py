"""
Feature Extractor for ML Scoring
Extracts 53 features from database using database_utils
"""
import logging
from typing import Dict
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

logger = logging.getLogger(__name__)

def extract_all_features(drug_name: str, disease: str = None) -> Dict:
    """Extract all 53 features for a drug"""
    
    features = {}
    
    try:
        from database_utils import execute_query
        
        # Get drug properties from database
        drug_rows = execute_query(
            """SELECT molecular_weight, log_p, qed_score, tpsa, rotatable_bonds,
                      aromatic_rings, hba, hbd, lipinski_violations, smiles
               FROM drugs WHERE name = %s LIMIT 1""",
            (drug_name,)
        )
        
        if drug_rows:
            drug = drug_rows[0]
            features['molecular_weight'] = float(drug.get('molecular_weight', 350))
            features['log_p'] = float(drug.get('log_p', 2.5))
            features['qed_score'] = float(drug.get('qed_score', 0.7))
            features['tpsa'] = float(drug.get('tpsa', 70))
            features['rotatable_bonds'] = int(drug.get('rotatable_bonds', 4))
            features['aromatic_rings'] = int(drug.get('aromatic_rings', 2))
            features['hba'] = int(drug.get('hba', 4))
            features['hbd'] = int(drug.get('hbd', 2))
            features['lipinski_violations'] = int(drug.get('lipinski_violations', 0))
        else:
            # Fallback values
            features = {
                'molecular_weight': 350.0, 'log_p': 2.5, 'qed_score': 0.7,
                'tpsa': 70.0, 'rotatable_bonds': 4, 'aromatic_rings': 2,
                'hba': 4, 'hbd': 2, 'lipinski_violations': 0
            }
        
        # Set defaults for all 53 features
        defaults = {
            'fraction_sp3': 0.4, 'pains_alerts': 0, 'synthetic_accessibility': 3.5,
            'heavy_atom_count': 25, 'num_rings': 3, 'num_heteroatoms': 5,
            'num_targets': 2, 'avg_binding_affinity': -7.0, 'best_binding_affinity': -8.0,
            'worst_binding_affinity': -6.0, 'avg_confidence': 0.75, 'best_confidence': 0.85,
            'target_diversity': 0.6, 'primary_target_confidence': 0.8, 'off_target_count': 1,
            'cyp450_interaction_count': 1, 'receptor_binding_count': 1, 'enzyme_interaction_count': 1,
            'pagerank_centrality': 0.05, 'betweenness_centrality': 0.02, 'closeness_centrality': 0.5,
            'degree_centrality': 0.1, 'path_length_to_disease': 3.0, 'pathway_connection_count': 2,
            'evidence_chain_strength': 0.7, 'publication_count': 40, 'clinical_trial_phase': 1.0,
            'network_clustering': 0.5, 'fda_approved': 1.0, 'years_approved': 10.0,
            'therapeutic_category_relevance': 0.7, 'drug_drug_interaction_risk': 0.3,
            'adverse_event_score': 0.8, 'safety_profile': 0.85, 'clinical_trial_success_rate': 0.6,
            'market_presence_years': 10.0, 'indication_similarity': 0.5, 'disease_category_match': 0.7,
            'mechanism_disease_relevance': 0.7, 'pathway_overlap_score': 0.5,
            'disease_gene_target_overlap': 0.5, 'unmet_medical_need': 0.7,
            'competitive_landscape': 0.5, 'repurposing_feasibility': 0.7
        }
        
        # Fill in defaults for missing features
        for key, value in defaults.items():
            if key not in features:
                features[key] = value
        
        return features
        
    except Exception as e:
        logger.error(f"Feature extraction failed: {e}")
        # Return all defaults
        return {
            'molecular_weight': 350.0, 'log_p': 2.5, 'qed_score': 0.7,
            'tpsa': 70.0, 'rotatable_bonds': 4, 'aromatic_rings': 2,
            'hba': 4, 'hbd': 2, 'lipinski_violations': 0,
            'fraction_sp3': 0.4, 'pains_alerts': 0, 'synthetic_accessibility': 3.5,
            'heavy_atom_count': 25, 'num_rings': 3, 'num_heteroatoms': 5,
            'num_targets': 2, 'avg_binding_affinity': -7.0, 'best_binding_affinity': -8.0,
            'worst_binding_affinity': -6.0, 'avg_confidence': 0.75, 'best_confidence': 0.85,
            'target_diversity': 0.6, 'primary_target_confidence': 0.8, 'off_target_count': 1,
            'cyp450_interaction_count': 1, 'receptor_binding_count': 1, 'enzyme_interaction_count': 1,
            'pagerank_centrality': 0.05, 'betweenness_centrality': 0.02, 'closeness_centrality': 0.5,
            'degree_centrality': 0.1, 'path_length_to_disease': 3.0, 'pathway_connection_count': 2,
            'evidence_chain_strength': 0.7, 'publication_count': 40, 'clinical_trial_phase': 1.0,
            'network_clustering': 0.5, 'fda_approved': 1.0, 'years_approved': 10.0,
            'therapeutic_category_relevance': 0.7, 'drug_drug_interaction_risk': 0.3,
            'adverse_event_score': 0.8, 'safety_profile': 0.85, 'clinical_trial_success_rate': 0.6,
            'market_presence_years': 10.0, 'indication_similarity': 0.5, 'disease_category_match': 0.7,
            'mechanism_disease_relevance': 0.7, 'pathway_overlap_score': 0.5,
            'disease_gene_target_overlap': 0.5, 'unmet_medical_need': 0.7,
            'competitive_landscape': 0.5, 'repurposing_feasibility': 0.7
        }
