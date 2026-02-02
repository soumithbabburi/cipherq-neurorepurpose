"""
Scoring Engine - IMPROVED WITH VALIDATION
More discriminating factors, disease-specific scoring
"""
import logging
import json

logger = logging.getLogger(__name__)

def score_drug(drug_name: str, disease: str = None, disease_pathway_count: int = 0) -> float:
    """
    Score drug with disease connection as PRIMARY factor
    
    Args:
        drug_name: Drug name
        disease: Target disease
        disease_pathway_count: Number of pathways connecting to disease (from pre-filter)
    
    Returns: 0.0 to 1.0 score
    """
    try:
        # Load data
        with open('drugs.json', 'r') as f:
            drugs = json.load(f)
        with open('drug_interactions.json', 'r') as f:
            interactions = json.load(f)
        
        # Find drug
        drug_lower = drug_name.lower()
        clean = drug_lower.replace(' hydrochloride', '').replace(' sodium', '').replace(', sterile', '')
        
        drug_key = None
        if drug_lower in drugs:
            drug_key = drug_lower
        elif clean in drugs:
            drug_key = clean
        else:
            for key in drugs.keys():
                if clean in key or key in clean:
                    drug_key = key
                    break
        
        if not drug_key:
            logger.warning(f"Drug {drug_name} not found")
            return 0.0
        
        drug_info = drugs[drug_key]
        score = 0.0
        
        # Component 1: DISEASE CONNECTION STRENGTH (0-0.40) - MOST IMPORTANT!
        if disease_pathway_count >= 5:
            score += 0.40  # Excellent disease connection
        elif disease_pathway_count >= 3:
            score += 0.30  # Good connection
        elif disease_pathway_count >= 1:
            score += 0.15  # Some connection
        # else: 0.0 - no connection, don't recommend!
        
        # Component 2: Target Confidence (0-0.30)
        if drug_key in interactions:
            targets = interactions[drug_key]
            sorted_targets = sorted(targets, key=lambda x: x.get('confidence_score', 0), reverse=True)
            
            if sorted_targets:
                top_conf = sorted_targets[0].get('confidence_score', 0)
                
                if top_conf >= 0.9:
                    score += 0.30
                elif top_conf >= 0.8:
                    score += 0.24
                elif top_conf >= 0.7:
                    score += 0.18
                elif top_conf >= 0.5:
                    score += 0.10
                else:
                    score += 0.05
        
        # Component 3: Molecular Quality (0-0.20)
        qed = drug_info.get('qed_score', 0)
        if qed >= 0.8:
            score += 0.20
        elif qed >= 0.6:
            score += 0.15
        elif qed >= 0.4:
            score += 0.10
        else:
            score += 0.05
        
        # Component 4: FDA Approval (0-0.10)
        if drug_info.get('approved'):
            score += 0.10
        
        logger.info(f"Score for {drug_name}: {score:.3f} (disease pathways: {disease_pathway_count})")
        return min(1.0, score)
        
    except Exception as e:
        logger.error(f"Scoring failed for {drug_name}: {e}")
        return 0.0


__all__ = ['score_drug']
