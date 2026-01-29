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
        """Get drugs by category based on their pathway profiles - NO HARDCODING!"""
        self._load_data()
        
        # Load pathways
        try:
            with open('pathways.json', 'r') as f:
                pathways_data = json.load(f)
            
            with open('protein_pathways.json', 'r') as f:
                protein_pathways = json.load(f)
        except:
            logger.warning("Could not load pathway files")
            pathways_data = {}
            protein_pathways = {}
        
        drugs = []
        category_lower = category.lower()
        
        # Category matching logic (based on pathway types)
        for drug_key, drug_info in self._drugs.items():
            # Must have interactions
            if drug_key not in self._interactions:
                continue
            
            # Get drug's targets
            targets = self._interactions[drug_key]
            
            # Get pathways for those targets
            drug_pathway_categories = set()
            for target in targets:
                gene = target.get('gene_symbol', '')
                if gene in protein_pathways:
                    pathway_ids = protein_pathways[gene]
                    for pw_id in pathway_ids[:10]:  # Check first 10 pathways
                        pw_info = pathways_data.get(pw_id, {})
                        pw_category = pw_info.get('category', '').lower()
                        pw_name = pw_info.get('name', '').lower()
                        
                        drug_pathway_categories.add(pw_category)
                        
                        # Also categorize based on pathway names
                        if any(x in pw_name for x in ['insulin', 'glucose', 'diabetes', 'metabolic']):
                            drug_pathway_categories.add('diabetic')
                        if any(x in pw_name for x in ['dopamin', 'parkinson']):
                            drug_pathway_categories.add('dopaminergic')
                        if any(x in pw_name for x in ['cardiac', 'cardiovascular', 'blood pressure']):
                            drug_pathway_categories.add('cardiovascular')
            
            # Check if drug matches requested category
            category_match = False
            
            if 'diabet' in category_lower:
                if 'diabetic' in drug_pathway_categories or 'metabolic' in drug_pathway_categories:
                    category_match = True
            elif 'cardiovascular' in category_lower or 'cardiac' in category_lower:
                if 'cardiovascular' in drug_pathway_categories:
                    category_match = True
            elif 'parkinson' in category_lower or 'dopamin' in category_lower:
                if 'dopaminergic' in drug_pathway_categories or 'neurological' in drug_pathway_categories:
                    category_match = True
            elif 'alzheimer' in category_lower:
                if 'neurological' in drug_pathway_categories:
                    category_match = True
            else:
                # Default: if it has interactions, include it
                category_match = True
            
            if category_match:
                drugs.append({
                    'name': drug_info.get('name', drug_key),
                    'category': ', '.join(drug_pathway_categories) if drug_pathway_categories else 'General',
                    'class': 'Small molecule',
                    'mechanism': f"Targets: {', '.join([t['gene_symbol'] for t in targets[:3]])}",
                    'smiles': drug_info.get('smiles'),
                    'qed_score': drug_info.get('qed_score', 0.5),
                    'fda_status': 'Approved' if drug_info.get('approved') else 'Unknown'
                })
            
            if len(drugs) >= limit:
                break
        
        logger.info(f"✅ Retrieved {len(drugs)} drugs for category '{category}'")
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
