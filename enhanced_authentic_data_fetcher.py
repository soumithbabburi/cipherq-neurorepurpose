import logging
import json
import pickle
import pandas as pd
import requests
from typing import List, Dict, Optional

logger = logging.getLogger(__name__)

# Load ML model globally for the scoring logic
try:
    with open('ml_scoring_model.pkl', 'rb') as f:
        ML_MODEL = pickle.load(f)
    logger.info("✅ ML scoring model loaded")
except Exception as e:
    logger.error(f"❌ Failed to load ML model: {e}")
    ML_MODEL = None

class EnhancedAuthenticDataFetcher:
    """
    Fetches real clinical trials and publications.
    Implements fallback strategies for metabolites and biological targets.
    """
    
    def __init__(self):
        self.clinicaltrials_base = "https://clinicaltrials.gov/api/v2/studies"
        self.pubmed_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        logger.info("Enhanced Data Fetcher initialized")

    def fetch_clinical_trials(self, drug_name: str, disease_name: str, max_results: int = 5) -> List[Dict]:
        """Fetch trials with fallback to metabolites if primary search fails."""
        trials = self._execute_trial_search(drug_name, disease_name, max_results)
        
        # Strategy: If no trials, check if we have metabolites in drugs.json to search for
        if not trials:
            try:
                with open('drugs.json', 'r') as f:
                    drugs = json.load(f)
                drug_info = drugs.get(drug_name.lower(), {})
                metabolites = drug_info.get('metabolites', [])
                for met in metabolites[:2]:
                    logger.info(f"Searching trials for metabolite: {met}")
                    trials = self._execute_trial_search(met, disease_name, max_results)
                    if trials: break
            except:
                pass
        return trials if trials else self._get_fallback_trials(drug_name, disease_name)

    def fetch_publications(self, drug_name: str, disease_name: str, max_results: int = 5) -> List[Dict]:
        """Fetch papers with fallback to targets (e.g., 'BACE1 inhibitor') if needed."""
        papers = self._execute_pubmed_search(f"{drug_name} AND {disease_name}", max_results)
        
        # Strategy: Target-Based Search if no direct publications found
        if not papers:
            try:
                with open('drug_interactions.json', 'r') as f:
                    interactions = json.load(f)
                targets = interactions.get(drug_name.lower(), [])
                if targets:
                    # Get highest confidence target
                    primary_target = max(targets, key=lambda x: x.get('confidence_score', 0))
                    gene = primary_target.get('gene_symbol', '')
                    if gene:
                        query = f"{gene} inhibitor AND {disease_name}"
                        logger.info(f"Searching target-based papers: {query}")
                        papers = self._execute_pubmed_search(query, max_results)
            except:
                pass
        return papers if papers else self._get_fallback_publications(drug_name, disease_name)

    def _execute_trial_search(self, drug, disease, max_results):
        """Internal helper for ClinicalTrials.gov API"""
        try:
            params = {"query.cond": disease, "query.term": drug, "pageSize": max_results, "format": "json"}
            response = requests.get(self.clinicaltrials_base, params=params, timeout=10)
            if response.status_code == 200:
                studies = response.json().get('studies', [])
                return [{
                    'nct_id': s.get('protocolSection', {}).get('identificationModule', {}).get('nctId'),
                    'title': s.get('protocolSection', {}).get('identificationModule', {}).get('briefTitle'),
                    'status': s.get('protocolSection', {}).get('statusModule', {}).get('overallStatus'),
                    'url': f"https://clinicaltrials.gov/study/{s.get('protocolSection', {}).get('identificationModule', {}).get('nctId')}"
                } for s in studies]
        except:
            return []
        return []

    def _execute_pubmed_search(self, query, max_results):
        """Internal helper for PubMed API"""
        try:
            search_params = {"db": "pubmed", "term": query, "retmax": max_results, "retmode": "json"}
            res = requests.get(f"{self.pubmed_base}/esearch.fcgi", params=search_params, timeout=10)
            ids = res.json().get('esearchresult', {}).get('idlist', [])
            if not ids: return []
            
            sum_res = requests.get(f"{self.pubmed_base}/esummary.fcgi", params={"db": "pubmed", "id": ','.join(ids), "retmode": "json"}, timeout=10)
            results = sum_res.json().get('result', {})
            return [{
                'pmid': pmid,
                'title': results[pmid].get('title'),
                'journal': results[pmid].get('source'),
                'year': results[pmid].get('pubdate', '')[:4],
                'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            } for pmid in ids if pmid in results]
        except:
            return []

    def _get_fallback_trials(self, drug, disease):
        return [{'nct_id': 'N/A', 'title': f'No specific trials found for {drug}', 'status': 'N/A', 'url': '#'}]

    def _get_fallback_publications(self, drug, disease):
        return [{'pmid': 'N/A', 'title': f'No specific publications for {drug}', 'url': '#'}]

def score_drug(drug_name: str, disease: str = None) -> float:
    """
    Strategy: In-Silico Prediction.
    Strictly uses the ML Model (100 decision trees) to score the drug.
    """
    if ML_MODEL is None:
        return 0.0
    try:
        with open('drugs.json', 'r') as f: drugs = json.load(f)
        with open('drug_interactions.json', 'r') as f: interactions = json.load(f)
        
        drug_key = drug_name.lower()
        if drug_key not in drugs: return 0.0
        
        # Features for ML: mw, target_conf, num_targets
        mw = drugs[drug_key].get('molecular_weight', 0)
        targets = interactions.get(drug_key, [])
        target_conf = max(t.get('confidence_score', 0) for t in targets) if targets else 0
        num_targets = len(targets)
        
        X = pd.DataFrame([[mw, target_conf, num_targets]], columns=['mw', 'target_conf', 'num_targets'])
        prob = ML_MODEL.predict_proba(X)[0][1]
        return float(prob)
    except:
        return 0.0
