"""
Scoring Engine - ML + Rule-based Hybrid
Uses ML model when available, falls back to rules
"""
import logging
import json

logger = logging.getLogger(__name__)

# Try to load ML model
ML_MODEL = None
try:
    import pickle
    with open('ml_scoring_model.pkl', 'rb') as f:
        ML_MODEL = pickle.load(f)
    logger.info("✅ ML scoring model loaded")
except:
    logger.info("ML model not available, using rule-based scoring")


def score_drug(drug_name: str, disease: str = None, disease_pathway_count: int = 0) -> float:
    """Score drug using ML (preferred) or rules (fallback)"""
    try:
        with open('drugs.json', 'r') as f:
            drugs = json.load(f)
        with open('drug_interactions.json', 'r') as f:
            interactions = json.load(f)
        
        drug_lower = drug_name.lower()
        clean = drug_lower.replace(' hydrochloride', '').replace(' sodium', '').replace(', sterile', '')
        
        drug_key = drug_lower if drug_lower in drugs else (clean if clean in drugs else None)
        
        if not drug_key:
            for key in drugs.keys():
                if clean in key or key in clean:
                    drug_key = key
                    break
        
        if not drug_key:
            return 0.0
        
        drug_info = drugs[drug_key]
        
        # TRY ML MODEL FIRST
        if ML_MODEL is not None:
            try:
                import pandas as pd
                
                mw = drug_info.get('molecular_weight', 0)
                
                if drug_key in interactions:
                    targets = interactions[drug_key]
                    target_conf = max(t.get('confidence_score', 0) for t in targets) if targets else 0
                    num_targets = len(targets)
                else:
                    target_conf = 0
                    num_targets = 0
                
                if mw > 0 and target_conf > 0:
                    X = pd.DataFrame([[mw, target_conf, num_targets]], 
                                   columns=['mw', 'target_conf', 'num_targets'])
                    prob = ML_MODEL.predict_proba(X)[0][1]
                    logger.info(f"ML score: {drug_name} = {prob:.3f}")
                    return float(prob)
            except Exception as ml_err:
                logger.warning(f"ML failed: {ml_err}")
        
        # RULE-BASED FALLBACK
        score = 0.0
        
        # Calculate pathways
        if disease_pathway_count == 0 and disease:
            try:
                with open('protein_pathways.json', 'r') as f:
                    protein_pathways = json.load(f)
                with open('pathways.json', 'r') as f:
                    pathways_data = json.load(f)
                
                if drug_key in interactions:
                    disease_lower = disease.lower()
                    for target in interactions[drug_key][:5]:
                        gene = target.get('gene_symbol', '').upper()
                        if gene in protein_pathways:
                            for pw_id in protein_pathways[gene]:
                                if pw_id in pathways_data:
                                    pw_name = pathways_data[pw_id].get('name', '').lower()
                                    if 'alzheimer' in disease_lower:
                                        if any(kw in pw_name for kw in ['insulin', 'glucose', 'metabol']):
                                            disease_pathway_count += 1
            except:
                pass
        
        # Pathway score
        if disease_pathway_count >= 5:
            score += 0.40
        elif disease_pathway_count >= 3:
            score += 0.30
        elif disease_pathway_count >= 1:
            score += 0.15
        
        # Target confidence
        if drug_key in interactions:
            targets = interactions[drug_key]
            if targets:
                top_conf = max(t.get('confidence_score', 0) for t in targets)
                if top_conf >= 0.9:
                    score += 0.30
                elif top_conf >= 0.7:
                    score += 0.20
                elif top_conf >= 0.5:
                    score += 0.10
        
        # QED
        qed = drug_info.get('qed_score', 0)
        if qed >= 0.6:
            score += 0.15
        elif qed >= 0.4:
            score += 0.10
        
        # FDA approval
        if drug_info.get('approved'):
            score += 0.10
        
        logger.info(f"Rule score: {drug_name} = {score:.3f}")
        return min(1.0, score)
        
    except Exception as e:
        logger.error(f"Scoring failed: {e}")
        return 0.0


__all__ = ['score_drug']
