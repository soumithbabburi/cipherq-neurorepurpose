"""
Workflow Optimizer - JSON ONLY VERSION
Uses scoring_engine.py (which now uses JSON)
"""
import logging
from typing import List, Dict
from scoring_engine import score_drug

logger = logging.getLogger(__name__)


def select_best_drugs_for_analysis(drugs: List, disease: str = None, limit: int = 3) -> List:
    """
    Select best drugs using JSON-based scoring
    
    Args:
        drugs: List of drug dicts or names
        disease: Target disease
        limit: Number of top drugs to return
    
    Returns:
        Top N drugs sorted by score from JSON data
    """
    try:
        # Score each drug
        scored_drugs = []
        
        for drug in drugs:
            drug_name = drug['name'] if isinstance(drug, dict) else str(drug)
            
            # Get score from JSON data
            ml_score = score_drug(drug_name, disease)
            
            # Add score to drug dict
            if isinstance(drug, dict):
                drug_copy = drug.copy()
                drug_copy['ml_score'] = ml_score
                drug_copy['confidence'] = ml_score
                scored_drugs.append(drug_copy)
            else:
                scored_drugs.append({
                    'name': drug_name,
                    'ml_score': ml_score,
                    'confidence': ml_score
                })
        
        # Sort by score (descending)
        scored_drugs.sort(key=lambda x: x.get('ml_score', 0), reverse=True)
        
        # Return top N
        top_drugs = scored_drugs[:limit]
        
        logger.info(f"✅ Selected top {len(top_drugs)} drugs by JSON-based score")
        if top_drugs:
            logger.info(f"   Top drug: {top_drugs[0]['name']} (score: {top_drugs[0]['ml_score']:.3f})")
        
        return top_drugs
        
    except Exception as e:
        logger.error(f"Drug selection failed: {e}")
        # Fallback to original order
        return drugs[:limit]


__all__ = ['select_best_drugs_for_analysis']
