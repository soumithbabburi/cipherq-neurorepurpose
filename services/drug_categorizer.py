import json
import os
import logging
from typing import List, Dict

# Relative import because they are in the same 'services' folder
from .drug_category_service import drug_category_service

logger = logging.getLogger(__name__)

class DrugCategorizer:
    def __init__(self):
        self._drugs = None
        self._interactions = None

    def _load_data(self):
        if self._drugs is not None:
            return
        try:
            # Get the folder where THIS python file is (app/services)
            current_dir = os.path.dirname(os.path.abspath(__file__))
            
            # Go UP one level to 'app' to find the JSONs
            drugs_path = os.path.normpath(os.path.join(current_dir, '..', 'drugs.json'))
            inter_path = os.path.normpath(os.path.join(current_dir, '..', 'drug_interactions.json'))

            with open(drugs_path, 'r') as f:
                self._drugs = json.load(f)
            with open(inter_path, 'r') as f:
                self._interactions = json.load(f)
            logger.info("✅ JSON Data loaded successfully")
        except Exception as e:
            logger.error(f"❌ Error loading JSON data from {current_dir}: {e}")
            self._drugs, self._interactions = {}, {}

    def get_drugs_by_category(self, category: str, limit: int = 50) -> List[Dict]:
        self._load_data()
        results = []
        target_cat = category.lower()

        for drug_key, drug_info in self._drugs.items():
            interactions = self._interactions.get(drug_key, [])
            genes = [i.get('gene_symbol') for i in interactions]

            # Use the service to categorize based on genes
            detected_cat, subcat = drug_category_service.categorize_by_genes(genes)

            if target_cat == 'general' or detected_cat.lower() == target_cat:
                results.append({
                    'name': drug_info.get('name', drug_key),
                    'category': detected_cat,
                    'subcategory': subcat,
                    'mechanism': interactions[0].get('protein_name', 'N/A') if interactions else 'N/A',
                    'smiles': drug_info.get('smiles'),
                    'fda_status': 'Approved' if drug_info.get('approved') else 'Unknown'
                })

            if len(results) >= limit:
                break
        return results

def get_drug_categorizer():
    return DrugCategorizer()
