"""
Disease Connection Filter
Pre-filters drugs to only include those with pathway connections to target disease
"""

import json
import logging

logger = logging.getLogger(__name__)


def filter_drugs_by_disease_connection(drug_list: list, disease_name: str) -> list:
    """
    Filter drugs to only include those with pathways connecting to target disease
    
    Args:
        drug_list: List of drug dicts with 'name' and 'targets'
        disease_name: Target disease (e.g., "Alzheimer's Disease")
    
    Returns:
        Filtered list with only drugs that have disease-relevant pathways
    """
    
    try:
        # Load pathway data
        with open('protein_pathways.json', 'r') as f:
            protein_pathways = json.load(f)
        
        with open('pathways.json', 'r') as f:
            pathways_data = json.load(f)
        
        disease_lower = disease_name.lower()
        filtered_drugs = []
        
        for drug in drug_list:
            drug_name = drug.get('name', '')
            targets = drug.get('targets', [])
            
            if not targets:
                continue
            
            # Check if ANY target has disease-relevant pathways
            has_disease_connection = False
            disease_pathway_count = 0
            
            for target_gene in targets[:5]:  # Check top 5 targets
                gene_upper = target_gene.upper()
                
                if gene_upper in protein_pathways:
                    pathway_ids = protein_pathways[gene_upper]
                    
                    # Check each pathway for disease relevance
                    for pw_id in pathway_ids:
                        if pw_id in pathways_data:
                            pathway_name = pathways_data[pw_id].get('name', '').lower()
                            
                            # Disease-specific matching
                            if 'alzheimer' in disease_lower:
                                if any(kw in pathway_name for kw in ['insulin', 'glucose', 'ampk', 'ppar', 'acetylcholine', 'cholin', 'metabol']):
                                    has_disease_connection = True
                                    disease_pathway_count += 1
                            
                            elif 'parkinson' in disease_lower:
                                if any(kw in pathway_name for kw in ['dopamin', 'mao', 'monoamine', 'motor']):
                                    has_disease_connection = True
                                    disease_pathway_count += 1
                            
                            elif 'diabetes' in disease_lower:
                                if any(kw in pathway_name for kw in ['insulin', 'glucose', 'metabol']):
                                    has_disease_connection = True
                                    disease_pathway_count += 1
                            
                            elif 'cardiovascular' in disease_lower or 'heart' in disease_lower:
                                if any(kw in pathway_name for kw in ['cardiac', 'vascular', 'blood pressure', 'lipid']):
                                    has_disease_connection = True
                                    disease_pathway_count += 1
            
            if has_disease_connection:
                # Add pathway count to drug info
                drug['disease_pathway_count'] = disease_pathway_count
                filtered_drugs.append(drug)
                logger.info(f"✓ {drug_name}: {disease_pathway_count} disease-relevant pathways")
            else:
                logger.info(f"✗ {drug_name}: No disease connections - filtered out")
        
        logger.info(f"FILTER RESULTS: {len(filtered_drugs)}/{len(drug_list)} drugs passed disease connection filter")
        
        return filtered_drugs
        
    except Exception as e:
        logger.error(f"Disease filter failed: {e}")
        # If filter fails, return original list (don't break the app)
        return drug_list


__all__ = ['filter_drugs_by_disease_connection']
