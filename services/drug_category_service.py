"""
Drug Category Service - Uses drug_therapeutic_categories.json
"""
import json
import logging

logger = logging.getLogger(__name__)

class DrugCategoryService:
    def __init__(self):
        self.categories_data = {}
        self._load_categories()
    
    def _load_categories(self):
        """Load pre-computed categories from JSON"""
        try:
            with open('drug_therapeutic_categories.json', 'r') as f:
                self.categories_data = json.load(f)
            logger.info(f"✅ Loaded categories for {len(self.categories_data)} drugs")
        except Exception as e:
            logger.error(f"❌ Error loading categories: {e}")
    
    def categorize_by_genes(self, gene_symbols):
        """
        Categorize based on gene symbols
        Returns: (main_category, subcategory)
        """
        gene_symbols = [g.upper() for g in gene_symbols]
        
        # Diabetic targets
        if any(g in ['PPARG', 'PPARA', 'DPP4', 'ABCC8', 'KCNJ11', 'SLC5A2', 'GLP1R', 'INSR'] for g in gene_symbols):
            return ('Metabolic', 'Diabetes/Insulin Regulator')
        
        # Cardiovascular
        if any(g in ['HMGCR', 'ACE', 'AGTR1', 'ADRB1', 'CACNA1C'] for g in gene_symbols):
            return ('Cardiovascular', 'Vascular Support')
        
        # Parkinson's
        if any(g in ['DRD2', 'DRD3', 'MAOB', 'COMT'] for g in gene_symbols):
            return ('Neurological', 'Parkinson\'s Disease')
        
        # Alzheimer's
        if any(g in ['ACHE', 'BCHE', 'GRIN1'] for g in gene_symbols):
            return ('Neurological', 'Alzheimer\'s Disease')
        
        return ('Other', 'Experimental')
    
    def categorize_all_drugs(self, drugs_list):
        """Categorize list of drugs"""
        # Not used in main workflow anymore
        return drugs_list
    
    def get_category_info(self, category):
        """Get category info"""
        return {'name': category, 'description': f'{category} drugs'}

# Singleton
drug_category_service = DrugCategoryService()
