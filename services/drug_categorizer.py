"""
Drug Categorizer - JSON ONLY VERSION
NO DATABASE - uses drugs.json and drug_interactions.json
"""
import logging
import json
from typing import List, Dict

logger = logging.getLogger(__name__)

class DrugCategorizer:
    """Categorize and retrieve drugs from JSON files"""
    
    def __init__(self):
        logger.info("Drug categorizer initialized (JSON only)")
        self._drugs = None
        self._interactions = None
    
    def _load_data(self):
        """Load JSON files once"""
        if self._drugs is not None:
            return
        
        try:
            with open('drugs.json', 'r') as f:
                self._drugs = json.load(f)
            logger.info(f"✅ Loaded {len(self._drugs)} drugs from JSON")
        except:
            self._drugs = {}
            logger.warning("drugs.json not found")
        
        try:
            with open('drug_interactions.json', 'r') as f:
                self._interactions = json.load(f)
            logger.info(f"✅ Loaded interactions for {len(self._interactions)} drugs")
        except:
            self._interactions = {}
            logger.warning("drug_interactions.json not found")
    
    def get_drugs_by_category(self, category: str, limit: int = 100) -> List[Dict]:
        """Get drugs by category from JSON"""
        self._load_data()
        
        drugs = []
        category_lower = category.lower()
        
        # Filter drugs that have interactions and match category
        for drug_key, drug_info in self._drugs.items():
            # Check if drug has targets (means it's useful)
            if drug_key not in self._interactions:
                continue
            
            # Simple category matching (can be improved)
            # For now, if drug has interactions, include it
            drugs.append({
                'name': drug_info.get('name', drug_key),
                'category': 'General',
                'class': 'Small molecule',
                'mechanism': 'Unknown',
                'smiles': drug_info.get('smiles'),
                'qed_score': drug_info.get('qed_score', 0.5),
                'fda_status': 'Approved' if drug_info.get('approved') else 'Unknown'
            })
            
            if len(drugs) >= limit:
                break
        
        logger.info(f"✅ Retrieved {len(drugs)} drugs")
        return drugs
    
    def get_all_categories(self) -> List[str]:
        """Get all categories"""
        # Return predefined categories since JSON doesn't have category field
        return [
            'Cardiovascular',
            'Metabolic', 
            'Anti-inflammatory',
            'Neuroprotective',
            'Psychiatric',
            'Antibiotic_Repurposed'
        ]
    
    def get_random_drugs(self, limit: int = 10) -> List[Dict]:
        """Get random drugs from JSON"""
        self._load_data()
        
        import random
        drug_items = list(self._drugs.items())
        random_items = random.sample(drug_items, min(limit, len(drug_items)))
        
        drugs = []
        for drug_key, drug_info in random_items:
            drugs.append({
                'name': drug_info.get('name', drug_key),
                'category': 'General',
                'class': 'Small molecule',
                'smiles': drug_info.get('smiles')
            })
        
        return drugs


# Singleton instance
_categorizer = None

def get_drug_categorizer():
    """Get singleton DrugCategorizer instance"""
    global _categorizer
    if _categorizer is None:
        _categorizer = DrugCategorizer()
    return _categorizer
