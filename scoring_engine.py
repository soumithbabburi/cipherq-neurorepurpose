import logging
import json
import pickle
import pandas as pd

logger = logging.getLogger(__name__)

# Load ML model - strict requirement
try:
    with open('ml_scoring_model.pkl', 'rb') as f:
        ML_MODEL = pickle.load(f)
    logger.info("✅ ML scoring model loaded")
except Exception as e:
    logger.error(f"❌ Failed to load ML model: {e}")
    ML_MODEL = None

def score_drug(drug_name: str) -> float:
    """Score drug using ML model only."""
    if ML_MODEL is None:
        logger.error("Scoring aborted: ML model not loaded.")
        return 0.0

    try:
        # Load necessary data
        with open('drugs.json', 'r') as f:
            drugs = json.load(f)
        with open('drug_interactions.json', 'r') as f:
            interactions = json.load(f)
        
        # Name normalization and lookup
        drug_lower = drug_name.lower()
        clean = drug_lower.replace(' hydrochloride', '').replace(' sodium', '').replace(', sterile', '')
        drug_key = drug_lower if drug_lower in drugs else (clean if clean in drugs else None)
        
        if not drug_key:
            for key in drugs.keys():
                if clean in key or key in clean:
                    drug_key = key
                    break
        
        if not drug_key:
            logger.warning(f"Drug '{drug_name}' not found in database.")
            return 0.0
        
        # Extract features for ML model [cite: 3, 207]
        drug_info = drugs[drug_key]
        mw = drug_info.get('molecular_weight', 0)
        
        targets = interactions.get(drug_key, [])
        target_conf = max(t.get('confidence_score', 0) for t in targets) if targets else 0
        num_targets = len(targets)
        
        # Ensure minimum feature requirements are met for valid prediction
        if mw > 0 and target_conf > 0:
            # Create DataFrame with exact feature names expected by the model [cite: 3, 207]
            X = pd.DataFrame([[mw, target_conf, num_targets]], 
                             columns=['mw', 'target_conf', 'num_targets'])
            
            # Predict probability (class 1)
            prob = ML_MODEL.predict_proba(X)[0][1]
            logger.info(f"ML score for {drug_name}: {prob:.3f}")
            return float(prob)
        
        logger.warning(f"Insufficient data for ML scoring: {drug_name}")
        return 0.0
            
    except Exception as e:
        logger.error(f"Scoring failed: {e}")
        return 0.0

__all__ = ['score_drug']
