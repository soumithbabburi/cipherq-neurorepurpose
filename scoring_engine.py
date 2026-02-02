"""
Scoring Engine - IMPROVED WITH VALIDATION
More discriminating factors, disease-specific scoring
"""
import logging
import json

logger = logging.getLogger(__name__)

def score_drug(drug_name: str, disease: str = None) -> float:
    """
    Score drug with improved discrimination
    Returns: 0.0 to 1.0 (more variation than before!)
    """
    try:
        # Load data
        with open('drugs.json', 'r') as f:
            drugs = json.load(f)
        with open('drug_interactions.json', 'r') as f:
            interactions = json.load(f)
        with open('protein_pathways.json', 'r') as f:
            protein_pathways = json.load(f)
        
        # Find drug with flexible matching
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
        
        # Component 1: QED score (0-0.25) - MORE DISCRIMINATING
        qed = drug_info.get('qed_score', 0)
        if qed >= 0.8:
            score += 0.25  # Excellent
        elif qed >= 0.6:
            score += 0.18  # Good
        elif qed >= 0.4:
            score += 0.10  # Fair
        else:
            score += 0.03  # Poor
        
        # Component 2: Target confidence (0-0.35) - MOST IMPORTANT
        if drug_key in interactions:
            targets = interactions[drug_key]
            
            # Sort by confidence
            sorted_targets = sorted(targets, key=lambda x: x.get('confidence_score', 0), reverse=True)
            
            # Top target confidence matters most
            if sorted_targets:
                top_conf = sorted_targets[0].get('confidence_score', 0)
                
                if top_conf >= 0.9:
                    score += 0.35  # Very high confidence
                elif top_conf >= 0.8:
                    score += 0.28  # High confidence
                elif top_conf >= 0.7:
                    score += 0.20  # Medium confidence
                elif top_conf >= 0.5:
                    score += 0.12  # Low confidence
                else:
                    score += 0.05  # Very low confidence
        
        # Component 3: Disease-specific pathways (0-0.25)
        if disease and drug_key in interactions:
            disease_lower = disease.lower()
            disease_relevant_pathways = 0
            
            for target in interactions[drug_key]:
                gene = target.get('gene_symbol', '')
                if gene in protein_pathways:
                    # Check if pathways are disease-relevant
                    for pathway_id in protein_pathways[gene]:
                        # Disease-specific pathway matching
                        if 'alzheimer' in disease_lower:
                            if any(kw in pathway_id.lower() for kw in ['insulin', 'glucose', 'metabolic', 'ampk', 'ppar', 'neuro', 'synap']):
                                disease_relevant_pathways += 1
                        elif 'parkinson' in disease_lower:
                            if any(kw in pathway_id.lower() for kw in ['dopamin', 'mao', 'catechol']):
                                disease_relevant_pathways += 1
                        elif 'diabetes' in disease_lower:
                            if any(kw in pathway_id.lower() for kw in ['insulin', 'glucose', 'metabolic']):
                                disease_relevant_pathways += 1
            
            if disease_relevant_pathways >= 5:
                score += 0.25
            elif disease_relevant_pathways >= 3:
                score += 0.18
            elif disease_relevant_pathways >= 1:
                score += 0.10
        
        # Component 4: FDA approval (0-0.15)
        if drug_info.get('approved'):
            score += 0.15
        
        logger.info(f"Score for {drug_name}: {score:.3f}")
        return min(1.0, score)
        
    except Exception as e:
        logger.error(f"Scoring failed for {drug_name}: {e}")
        return 0.0


__all__ = ['score_drug']
