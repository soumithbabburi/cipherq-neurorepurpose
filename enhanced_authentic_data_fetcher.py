import logging
import json
import pickle
import pandas as pd
import requests
from typing import List, Dict, Optional

logger = logging.getLogger(__name__)

# Load ML model
try:
    with open('ml_scoring_model.pkl', 'rb') as f:
        ML_MODEL = pickle.load(f)
    logger.info("✅ ML scoring model loaded")
except Exception as e:
    logger.error(f"❌ Failed to load ML model: {e}")
    ML_MODEL = None

class SmartDataFetcher:
    """Enhanced fetcher with fallback strategies for metabolites and targets."""
    
    def __init__(self):
        self.clinicaltrials_base = "https://clinicaltrials.gov/api/v2/studies"
        self.pubmed_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def fetch_comprehensive_evidence(self, drug_name: str, disease_name: str, drug_info: Dict, targets: List[Dict]) -> Dict:
        """
        Implements multi-stage search:
        1. Direct (Drug + Disease)
        2. Analogs/Metabolites (from drugs.json)
        3. Target-based (e.g., 'BACE1 inhibitors')
        """
        # Strategy 1: Direct Search
        trials = self._fetch_trials(drug_name, disease_name)
        papers = self._fetch_papers(drug_name, disease_name)
        
        # Strategy 2: Metabolites/Analogs Fallback
        if not papers or not trials:
            logger.info(f"Low evidence for {drug_name}. Searching for metabolites/synonyms...")
            alternatives = drug_info.get('metabolites', []) + drug_info.get('synonyms', [])
            for alt in alternatives[:2]: # Check top 2 alternatives
                if not trials: trials = self._fetch_trials(alt, disease_name)
                if not papers: papers = self._fetch_papers(alt, disease_name)

        # Strategy 3: Target-Based Search (e.g., 'MAPK1 inhibitor')
        if not papers and targets:
            primary_target = max(targets, key=lambda x: x.get('confidence_score', 0))
            gene = primary_target.get('gene_symbol', '')
            if gene:
                target_query = f"{gene} inhibitor"
                logger.info(f"No direct papers. Searching target: {target_query}")
                papers = self._fetch_papers(target_query, disease_name)

        return {
            "trials": trials,
            "publications": papers,
            "has_real_world_data": len(trials) > 0 or len(papers) > 0
        }

    def _fetch_trials(self, query_term: str, disease: str) -> List:
        # (Standard ClinicalTrials.gov API request logic here...)
        return [] # Placeholder for actual API implementation

    def _fetch_papers(self, query_term: str, disease: str) -> List:
        # (Standard PubMed API request logic here...)
        return [] # Placeholder for actual API implementation


def score_drug(drug_name: str, disease: str) -> Dict:
    """
    Final Scoring Engine:
    - Uses ML Model as the ground truth.
    - Uses SmartDataFetcher for supporting evidence.
    - Strategy 3: Relies solely on ML if Real-World Data (RWD) is missing.
    """
    try:
        # 1. Load Data
        with open('drugs.json', 'r') as f: drugs = json.load(f)
        with open('drug_interactions.json', 'r') as f: interactions = json.load(f)
        
        drug_key = drug_name.lower() # Simplified lookup
        if drug_key not in drugs: return {"error": "Drug not found"}
        
        drug_info = drugs[drug_key]
        drug_targets = interactions.get(drug_key, [])
        
        # 2. ML FEATURE EXTRACTION
        mw = drug_info.get('molecular_weight', 0)
        target_conf = max([t.get('confidence_score', 0) for t in drug_targets]) if drug_targets else 0
        num_targets = len(drug_targets)
        
        # 3. CALCULATE ML SCORE (In-Silico Prediction)
        ml_score = 0.0
        if ML_MODEL and mw > 0:
            X = pd.DataFrame([[mw, target_conf, num_targets]], 
                             columns=['mw', 'target_conf', 'num_targets'])
            ml_score = float(ML_MODEL.predict_proba(X)[0][1])

        # 4. FETCH EVIDENCE (Metabolite & Target strategies included)
        fetcher = SmartDataFetcher()
        evidence = fetcher.fetch_comprehensive_evidence(drug_name, disease, drug_info, drug_targets)
        
        # 5. INTEGRATION LOGIC
        # Strategy: If no real-world evidence exists, the ML score IS the final score.
        final_score = ml_score
        status = "Predicted (In-Silico)" if not evidence['has_real_world_data'] else "Validated (RWD)"

        return {
            "drug": drug_name,
            "score": round(final_score, 4),
            "confidence_type": status,
            "evidence": evidence
        }

    except Exception as e:
        logger.error(f"Scoring failed: {e}")
        return {"score": 0.0, "error": str(e)}
