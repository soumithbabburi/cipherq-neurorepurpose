import json
import os
import logging
from typing import List, Dict

# Relative import
from .drug_category_service import drug_category_service

logger = logging.getLogger(__name__)

class DrugCategorizer:
    def __init__(self):
        self._drugs = None
        self._interactions = None
        
        # MAP UI CATEGORIES TO INTERNAL LOGIC
        # If user clicks "Diabetic Drugs", we search for "Metabolic"
        self.category_mapping = {
            'diabetic drugs': 'metabolic',
            'diabetes': 'metabolic',
            'cardiovascular': 'cardiovascular',
            'neuroprotective': 'neuroprotective',
            'anti-inflammatory': 'anti-inflammatory',
            'psychiatric': 'psychiatric'
        }

    def _load_data(self):
        if self._drugs is not None:
            return
        try:
            # Go up one level to find JSONs
            current_dir = os.path.dirname(os.path.abspath(__file__))
            drugs_path = os.path.normpath(os.path.join(current_dir, '..', 'drugs.json'))
            inter_path = os.path.normpath(os.path.join(current_dir, '..', 'drug_interactions.json'))

            with open(drugs_path, 'r') as f:
                self._drugs = json.load(f)
            with open(inter_path, 'r') as f:
                self._interactions = json.load(f)
        except Exception as e:
            logger.error(f"❌ Error loading JSONs: {e}")
            self._drugs, self._interactions = {}, {}

    def get_drugs_by_category(self, category: str, limit: int = 50) -> List[Dict]:
        self._load_data()
        results = []
        
        # 1. Normalize Category Name
        # If user asks for "Diabetic Drugs", convert it to "metabolic"
        query_cat = category.lower()
        target_cat = self.category_mapping.get(query_cat, query_cat)

        logger.info(f"🔍 Searching for category: '{query_cat}' -> mapped to '{target_cat}'")

        for drug_key, drug_info in self._drugs.items():
            interactions = self._interactions.get(drug_key, [])
            genes = [i.get('gene_symbol') for i in interactions]

            # 2. Get Category from Service
            detected_cat, subcat = drug_category_service.categorize_by_genes(genes)
            
            # 3. Match Logic
            # matches "Metabolic" == "Metabolic"
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
                
        logger.info(f"✅ Found {len(results)} drugs for '{category}'")
        return results

def get_drug_categorizer():
    return DrugCategorizer()
