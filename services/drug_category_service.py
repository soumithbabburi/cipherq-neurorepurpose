import pandas as pd
import logging
from typing import List, Tuple

logger = logging.getLogger(__name__)

class DrugCategoryService:
    def __init__(self):
        self.therapeutic_taxonomy = {
            'Cardiovascular': {'relevance_score': 0.72},
            'Metabolic': {'relevance_score': 0.68},
            'Anti-inflammatory': {'relevance_score': 0.65},
            'Neuroprotective': {'relevance_score': 0.78},
            'Psychiatric': {'relevance_score': 0.61}
        }
        self.gene_metadata = {}
        self._load_tsv_metadata()

    def _load_tsv_metadata(self):
        """Processes the categories (1).tsv file to map genes to functions"""
        try:
            # Note: Ensure this file is uploaded to your GitHub repo
            df = pd.read_csv('categories (1).tsv', sep='\t')
            self.gene_metadata = df.groupby('name')['name-2'].apply(set).to_dict()
            logger.info("✅ TSV Metadata loaded successfully")
        except Exception as e:
            logger.error(f"❌ Error loading TSV: {e}")

    def categorize_by_genes(self, gene_symbols: List[str]) -> Tuple[str, str]:
        """Determines the category (e.g., Metabolic) by checking target genes"""
        functions = set()
        for gene in gene_symbols:
            if gene in self.gene_metadata:
                functions.update(self.gene_metadata[gene])

        # DIABETES / METABOLIC RULES
        metabolic_markers = {'DPP4', 'PPARG', 'ABCC8', 'KCNJ11', 'SLC5A2', 'PRKAA1', 'PRKAA2', 'GLP1R'}
        if any(g in metabolic_markers for g in gene_symbols):
            return 'Metabolic', 'Insulin & Glucose Regulation'

        # CARDIOVASCULAR RULES
        if 'ION CHANNEL' in functions or any(g.startswith(('ADR', 'ACE')) for g in gene_symbols):
            return 'Cardiovascular', 'Vascular Support'

        # ANTI-INFLAMMATORY RULES
        if 'PROTEASE' in functions or any(g in {'PTGS2', 'TNF', 'IL1B'} for g in gene_symbols):
            return 'Anti-inflammatory', 'Neuroinflammation Control'

        return 'Other', 'Experimental'

# Create the global instance that drug_categorizer.py will import
drug_category_service = DrugCategoryService()
