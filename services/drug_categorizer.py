import json
import logging
from typing import List, Dict
# This line looks for drug_category_service.py in the same folder
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
        except Exception as e:
            logger.error(f"Error loading JSON data: {e}")
            self._drugs, self._interactions = {}, {}

    def get_drugs_by_category(self, category: str, limit: int = 50) -> List[Dict]:
        self._load_data()
        results = []
        target_cat = category.lower()

        for drug_key, drug_info in self._drugs.items():
            interactions = self._interactions.get(drug_key, [])
            genes = [i.get('gene_symbol') for i in interactions]

            # Use the service to find the category based on the TSV + Genes
            detected_cat, subcat = drug_category_service.categorize_by_genes(genes)

            if target_cat == 'general' or detected_cat.lower() == target_cat:
                results.append({
                    'name': drug_info.get('name', drug_key),
                    'category': detected_cat,
                    'subcategory': subcat,
                    'smiles': drug_info.get('smiles'),
                    'fda_status': 'Approved' if drug_info.get('approved') else 'Unknown'
                })

            if len(results) >= limit:
                break
        return results

# Interface for your main app
def get_drug_categorizer():
    return DrugCategorizer()
