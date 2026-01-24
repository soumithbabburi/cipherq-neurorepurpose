"""
Workflow Optimizer - ML-powered drug selection
"""
import logging
from typing import List, Dict
from scoring_engine import score_drug

logger = logging.getLogger(__name__)


def select_best_drugs_for_analysis(drugs: List, disease: str = None, limit: int = 3) -> List:
    """
    Select best drugs using ML scoring
    
    Args:
        drugs: List of drug dicts or names
        disease: Target disease
        limit: Number of top drugs to return
    
    Returns:
        Top N drugs sorted by ML score
    """
    try:
        # Score each drug
        scored_drugs = []
        
        for drug in drugs:
            drug_name = drug['name'] if isinstance(drug, dict) else str(drug)
            
            # Get ML score
            ml_score = score_drug(drug_name, disease)
            
            # Add score to drug dict
            if isinstance(drug, dict):
                drug_copy = drug.copy()
                drug_copy['ml_score'] = ml_score
                scored_drugs.append(drug_copy)
            else:
                scored_drugs.append({
                    'name': drug_name,
                    'ml_score': ml_score
                })
        
        # Sort by ML score (descending)
        scored_drugs.sort(key=lambda x: x.get('ml_score', 0), reverse=True)
        
        # Return top N
        top_drugs = scored_drugs[:limit]
        
        logger.info(f"✅ Selected top {len(top_drugs)} drugs by ML score")
        
        return top_drugs
        
    except Exception as e:
        logger.error(f"Drug selection failed: {e}")
        # Fallback to confidence-based sorting
        if isinstance(drugs[0], dict):
            sorted_drugs = sorted(drugs, key=lambda x: x.get('confidence', 0.5), reverse=True)
        else:
            sorted_drugs = drugs
        
        return sorted_drugs[:limit]


__all__ = ['select_best_drugs_for_analysis']