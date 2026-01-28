import json
import logging
from typing import List, Dict
from drug_category_service import drug_category_service

logger = logging.getLogger(__name__)

class DrugCategorizer:
    def __init__(self):
        self._drugs = None
        self._interactions = None

    def _load_data(self):
        if self._drugs is not None:
            return
        
        try:
            with open('drugs.json', 'r') as f:
                self._drugs = json.load(f)
            with open('drug_interactions.json', 'r') as f:
                self._interactions = json.load(f)
            logger.info("JSON data loaded into Categorizer")
        except Exception as e:
            logger.error(f"Error loading JSONs: {e}")
            self._drugs, self._interactions = {}, {}

    def get_drugs_by_category(self, category: str, limit: int = 100) -> List[Dict]:
        self._load_data()
        results = []
        category_lower = category.lower()

        for drug_key, drug_info in self._drugs.items():
            # Get interaction genes for this drug
            interactions = self._interactions.get(drug_key, [])
            gene_symbols = [i.get('gene_symbol') for i in interactions]

            # Use the TSV-based service to get the real category
            detected_cat, subcat = drug_category_service.categorize_by_genes(gene_symbols)

            # Filter logic
            if category_lower == 'general' or detected_cat.lower() == category_lower:
                results.append({
                    'name': drug_info.get('name', drug_key),
                    'category': detected_cat,
                    'subcategory': subcat,
                    'mechanism': interactions[0].get('protein_name', 'Unknown') if interactions else 'Unknown',
                    'smiles': drug_info.get('smiles'),
                    'qed_score': drug_info.get('qed_score', 0.5),
                    'fda_status': 'Approved' if drug_info.get('approved') else 'Unknown'
                })

            if len(results) >= limit:
                break
        
        return results

# Singleton instance
_categorizer = DrugCategorizer()

def get_drug_categorizer():
    return _categorizer
