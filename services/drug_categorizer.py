"""
Drug Categorizer - JSON ONLY VERSION
Corrected to filter drugs by category using keyword matching.
"""
import logging
import json
import random
from typing import List, Dict

logger = logging.getLogger(__name__)

class DrugCategorizer:
    """Categorize and retrieve drugs from JSON files"""
    
    def __init__(self):
        logger.info("Drug categorizer initialized (JSON only)")
        self._drugs = None
        self._interactions = None
        
        # Define keyword mappings for categories
        # This allows us to "guess" the category since the JSON lacks a category field.
        self.category_keywords = {
            'cardiovascular': ['ace', 'beta', 'blocker', 'statin', 'heart', 'blood pressure', 'angiotensin', 'adrenoceptor'],
            'metabolic': ['insulin', 'diabetes', 'glucose', 'metformin', 'lipid', 'glucagon'],
            'anti-inflammatory': ['cox', 'prostaglandin', 'aspirin', 'steroid', 'ibuprofen', 'interleukin', 'tnf'],
            'neuroprotective': ['glutamate', 'oxidative', 'nerve', 'amyloid'],
            'psychiatric': ['serotonin', 'dopamine', 'gaba', 'antipsychotic', 'depressant', 'reuptake'],
            'antibiotic_repurposed': ['bacterial', 'penicillin', 'enzyme', 'ribosome', 'antibiotic']
        }
    
    def _load_data(self):
        """Load JSON files once"""
        if self._drugs is not None:
            return
        
        try:
            with open('drugs.json', 'r') as f:
                self._drugs = json.load(f)
            logger.info(f"✅ Loaded {len(self._drugs)} drugs from JSON")
        except FileNotFoundError:
            self._drugs = {}
            logger.warning("drugs.json not found")
        
        try:
            with open('drug_interactions.json', 'r') as f:
                self._interactions = json.load(f)
            logger.info(f"✅ Loaded interactions for {len(self._interactions)} drugs")
        except FileNotFoundError:
            self._interactions = {}
            logger.warning("drug_interactions.json not found")
    
    def get_drugs_by_category(self, category: str, limit: int = 100) -> List[Dict]:
        """Get drugs filtered by category logic"""
        self._load_data()
        
        drugs = []
        category_lower = category.lower()
        keywords = self.category_keywords.get(category_lower, [])
        
        for drug_key, drug_info in self._drugs.items():
            # 1. Skip if no interaction data exists (ensures we have biological context)
            if drug_key not in self._interactions:
                continue
            
            # 2. Get target data for this drug
            drug_interactions = self._interactions[drug_key]
            targets = [inter.get('protein_name', '').lower() for inter in drug_interactions]
            drug_name = drug_info.get('name', drug_key).lower()
            
            # 3. Determine if drug belongs in the category
            # We match keywords against the drug name or its protein targets
            is_match = any(kw in drug_name or any(kw in t for t in targets) for kw in keywords)
            
            # If "General" is chosen, or if it's a specific match, add it
            if category_lower == 'general' or is_match:
                drugs.append({
                    'name': drug_info.get('name', drug_key),
                    'category': category.capitalize(),
                    'class': 'Small molecule',
                    'mechanism': targets[0].capitalize() if targets else 'Unknown',
                    'smiles': drug_info.get('smiles'),
                    'qed_score': drug_info.get('qed_score', 0.5),
                    'fda_status': 'Approved' if drug_info.get('approved') else 'Unknown'
                })
            
            if len(drugs) >= limit:
                break
        
        logger.info(f"✅ Retrieved {len(drugs)} drugs for category: {category}")
        return drugs
    
    def get_all_categories(self) -> List[str]:
        """Get all available categories"""
        return ['General'] + [k.capitalize() for k in self.category_keywords.keys()]
    
    def get_random_drugs(self, limit: int = 10) -> List[Dict]:
        """Get random drugs from JSON"""
        self._load_data()
        
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
