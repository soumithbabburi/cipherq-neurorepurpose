"""
Scoring Engine - JSON ONLY VERSION
NO DATABASE - uses drugs.json, drug_interactions.json, protein_pathways.json
"""
import logging
import json

logger = logging.getLogger(__name__)

def score_drug(drug_name: str, disease: str = None) -> float:
    """
    Score a drug using JSON data
    Returns: 0.0 to 1.0 confidence score
    """
    try:
        # Load drugs.json
        with open('drugs.json', 'r') as f:
            drugs = json.load(f)
        
        # Load drug_interactions.json
        with open('drug_interactions.json', 'r') as f:
            interactions = json.load(f)
        
        # Load protein_pathways.json
        with open('protein_pathways.json', 'r') as f:
            protein_pathways = json.load(f)
        
        # Try to find drug (flexible matching)
        drug_lower = drug_name.lower()
        drug_key = None
        
        # Exact match
        if drug_lower in drugs:
            drug_key = drug_lower
        else:
            # Partial match
            for key in drugs.keys():
                if drug_lower in key or key in drug_lower:
                    drug_key = key
                    break
        
        if not drug_key:
            logger.warning(f"Drug {drug_name} not found in JSON")
            return 0.0
        
        drug_info = drugs[drug_key]
        
        # Calculate score (0-1)
        score = 0.0
        
        # Component 1: QED score (0-0.3)
        qed = drug_info.get('qed_score', 0)
        if qed:
            score += float(qed) * 0.3
        
        # Component 2: Has targets (0-0.3)
        if drug_key in interactions:
            targets = interactions[drug_key]
            high_conf = sum(1 for t in targets if t.get('confidence_score', 0) > 0.8)
            score += min(0.3, high_conf * 0.1)
        
        # Component 3: Has pathways (0-0.2)
        if drug_key in interactions:
            pathway_count = 0
            for target in interactions[drug_key]:
                gene = target.get('gene_symbol', '')
                if gene in protein_pathways:
                    pathway_count += len(protein_pathways[gene])
            
            if pathway_count > 10:
                score += 0.2
            elif pathway_count > 5:
                score += 0.15
            elif pathway_count > 0:
                score += 0.1
        
        # Component 4: FDA approved (0-0.2)
        if drug_info.get('approved'):
            score += 0.2
        
        logger.info(f"Score for {drug_name}: {score:.3f}")
        return round(score, 3)
        
    except Exception as e:
        logger.error(f"Scoring failed for {drug_name}: {e}")
        return 0.0


def rank_drugs_for_disease(drugs: list, disease: str) -> list:
    """
    Rank drugs by ML score using JSON data
    
    Args:
        drugs: List of drug dicts or names
        disease: Target disease
    
    Returns:
        Sorted list of drugs
    """
    try:
        scored_drugs = []
        
        for drug in drugs:
            drug_name = drug['name'] if isinstance(drug, dict) else str(drug)
            
            # Calculate score
            ml_score = score_drug(drug_name, disease)
            
            if isinstance(drug, dict):
                drug_copy = drug.copy()
                drug_copy['confidence'] = ml_score
                drug_copy['ml_score'] = ml_score
                scored_drugs.append(drug_copy)
            else:
                scored_drugs.append({
                    'name': drug_name,
                    'confidence': ml_score,
                    'ml_score': ml_score
                })
        
        # Sort by score
        scored_drugs.sort(key=lambda x: x.get('ml_score', 0), reverse=True)
        
        logger.info(f"✅ Ranked {len(scored_drugs)} drugs for {disease}")
        
        return scored_drugs
        
    except Exception as e:
        logger.error(f"Ranking failed: {e}")
        return drugs


__all__ = ['score_drug', 'rank_drugs_for_disease']
