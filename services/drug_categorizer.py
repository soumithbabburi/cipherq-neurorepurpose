"""
Drug Categorizer
Uses drug_therapeutic_categories.json for proper categorization
"""
import logging
import json
from typing import List, Dict

logger = logging.getLogger(__name__)

class DrugCategorizer:
    """Categorize and retrieve drugs using pre-computed categories"""
    
    def __init__(self):
        logger.info("Drug categorizer initialized (using drug_therapeutic_categories.json)")
        self._drugs = None
        self._interactions = None
        self._categories = None
    
    def _load_data(self):
        """Load JSON files once"""
        if self._drugs is not None:
            return
        
        try:
            with open('drugs.json', 'r') as f:
                self._drugs = json.load(f)
            logger.info(f"✅ Loaded {len(self._drugs)} drugs")
        except:
            self._drugs = {}
        
        try:
            with open('drug_interactions.json', 'r') as f:
                self._interactions = json.load(f)
            logger.info(f"✅ Loaded interactions for {len(self._interactions)} drugs")
        except:
            self._interactions = {}
        
        try:
            with open('drug_therapeutic_categories.json', 'r') as f:
                self._categories = json.load(f)
            logger.info(f"✅ Loaded categories for {len(self._categories)} drugs")
        except:
            self._categories = {}
            logger.warning("drug_therapeutic_categories.json not found!")
    
    def get_drugs_by_category(self, category: str, limit: int = 100) -> List[Dict]:
        """Get drugs by therapeutic category"""
        self._load_data()
        
        if not self._categories:
            logger.error("Categories not loaded!")
            return []
        
        drugs = []
        
        # Normalize category
        category_map = {
            'diabetes': 'Diabetic',
            'diabetic': 'Diabetic',
            'cardiovascular': 'Cardiovascular',
            'parkinson': 'Parkinsons',
            'alzheimer': 'Alzheimers',
            'cancer': 'Cancer',
            'pain': 'Pain/Anti-inflammatory'
        }
        
        target_category = category_map.get(category.lower(), category)
        
        # Filter drugs
        for drug_name, drug_cats in self._categories.items():
            if target_category in drug_cats:
                drug_info = self._drugs.get(drug_name, {})
                targets_info = self._interactions.get(drug_name, [])
                
                drugs.append({
                    'name': drug_info.get('name', drug_name),
                    'category': ', '.join(drug_cats),
                    'class': 'Small molecule',
                    'mechanism': f"Targets: {', '.join([t['gene_symbol'] for t in targets_info[:3]])}",
                    'smiles': drug_info.get('smiles'),
                    'qed_score': drug_info.get('qed_score', 0.5),
                    'fda_status': 'Approved' if drug_info.get('approved') else 'Unknown',
                    'targets': [t['gene_symbol'] for t in targets_info]
                })
            
            if len(drugs) >= limit:
                break
        
        logger.info(f"✅ Retrieved {len(drugs)} drugs for '{category}'")
        return drugs
    
    def get_all_categories(self) -> List[str]:
        """Get all categories"""
        return [
            'Diabetic', 'Cardiovascular', 'Parkinsons', 'Alzheimers',
            'Cancer', 'Pain/Anti-inflammatory', 'Psychiatric', 'General'
        ]
    
    def get_random_drugs(self, limit: int = 10) -> List[Dict]:
        """Get random drugs"""
        self._load_data()
        
        import random
        drug_items = list(self._drugs.items())
        random_items = random.sample(drug_items, min(limit, len(drug_items)))
        
        return [{'name': d[1].get('name', d[0]), 'smiles': d[1].get('smiles')} for d in random_items]


# Singleton
_categorizer = None

def get_drug_categorizer():
    global _categorizer
    if _categorizer is None:
        _categorizer = DrugCategorizer()
    return _categorizer
