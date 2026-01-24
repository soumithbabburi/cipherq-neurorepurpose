"""
Scoring Engine - Simple Wrapper
Uses ML ensemble under the hood
"""
import logging

logger = logging.getLogger(__name__)

# Try to use ML scoring engine
try:
    from ml_scoring_engine_ENHANCED import score_drug as ml_score_drug
    from ml_scoring_engine_ENHANCED import rank_drugs_for_disease as ml_rank_drugs
    from ml_scoring_engine_ENHANCED import ml_scoring_engine
    
    ML_AVAILABLE = True
    logger.info("✅ ML Scoring Engine loaded (LightGBM + XGBoost + RF ensemble)")
    
    def score_drug(drug_name: str, disease: str = None) -> float:
        """Score a drug using ML ensemble"""
        return ml_score_drug(drug_name, disease)
    
    def rank_drugs_for_disease(drugs: list, disease: str) -> list:
        """Rank drugs using ML ensemble"""
        return ml_rank_drugs(drugs, disease)
    
except ImportError as e:
    logger.warning(f"ML scoring not available: {e}")
    logger.warning("Using database QED scores instead")
    ML_AVAILABLE = False
    
    def score_drug(drug_name: str, disease: str = None) -> float:
        """Get score from database QED - NO HARDCODING"""
        try:
            import psycopg2
            import os
            
            conn = psycopg2.connect(
                host=os.getenv("DB_HOST", "localhost"),
                database=os.getenv("DB_NAME", "cipherq_repurpose"),
                user=os.getenv("DB_USER", "babburisoumith"),
                password=os.getenv("DB_PASSWORD", "")
            )
            cursor = conn.cursor()
            cursor.execute("SELECT qed_score FROM drugs WHERE name = %s LIMIT 1", (drug_name,))
            result = cursor.fetchone()
            cursor.close()
            conn.close()
            
            if result and result[0]:
                logger.info(f"Score for {drug_name}: {result[0]} (from database)")
                return float(result[0])
            else:
                logger.warning(f"No score in database for {drug_name}")
                return 0.0  # Return 0 to indicate no score, not fake score
        except Exception as e:
            logger.error(f"Score lookup failed: {e}")
            return 0.0
    
    def rank_drugs_for_disease(drugs: list, disease: str) -> list:
        """Rank by database scores - NO FAKE SCORES"""
        for drug in drugs:
            if 'confidence' not in drug or not drug['confidence']:
                score = score_drug(drug.get('name', ''), disease)
                drug['confidence'] = score
        
        return sorted(drugs, key=lambda x: x.get('confidence', 0.0), reverse=True)


__all__ = ['score_drug', 'rank_drugs_for_disease']