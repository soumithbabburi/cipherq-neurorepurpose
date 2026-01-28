import pandas as pd
import os
import logging
from typing import List, Tuple

logger = logging.getLogger(__name__)

class DrugCategoryService:
    def __init__(self):
        self.gene_metadata = {}
        self._load_tsv_metadata()

    def _load_tsv_metadata(self):
        """Processes the categories (1).tsv file located in the parent folder"""
        try:
            # Get the folder where THIS python file is (app/services)
            current_dir = os.path.dirname(os.path.abspath(__file__))
            
            # Go UP one level to 'app' to find the TSV
            # path becomes: .../app/categories (1).tsv
            tsv_path = os.path.join(current_dir, '..', 'categories (1).tsv')
            
            # Normalize path (fixes slash issues)
            tsv_path = os.path.normpath(tsv_path)
            
            if not os.path.exists(tsv_path):
                logger.error(f"❌ TSV file not found at: {tsv_path}")
                return

            df = pd.read_csv(tsv_path, sep='\t')
            # Create a dictionary for gene function lookup
            self.gene_metadata = df.groupby('name')['name-2'].apply(set).to_dict()
            logger.info(f"✅ TSV Metadata loaded successfully from: {tsv_path}")
        except Exception as e:
            logger.error(f"❌ Error loading TSV: {e}")

    def categorize_by_genes(self, gene_symbols: List[str]) -> Tuple[str, str]:
        """Determines category based on TSV and Gene rules"""
        functions = set()
        for gene in gene_symbols:
            if gene in self.gene_metadata:
                functions.update(self.gene_metadata[gene])

        # METABOLIC RULES (Diabetes)
        metabolic_markers = {'DPP4', 'PPARG', 'ABCC8', 'KCNJ11', 'SLC5A2', 'PRKAA1', 'PRKAA2', 'GLP1R'}
        if any(g in metabolic_markers for g in gene_symbols):
            return 'Metabolic', 'Insulin & Glucose Regulation'

        # CARDIOVASCULAR RULES
        if 'ION CHANNEL' in functions or any(g.startswith(('ADR', 'ACE')) for g in gene_symbols):
            return 'Cardiovascular', 'Vascular Support'

        # ANTI-INFLAMMATORY RULES
        if 'PROTEASE' in functions or any(g in {'PTGS2', 'TNF', 'IL1B'} for g in gene_symbols):
            return 'Anti-inflammatory', 'Neuroinflammation Control'

        # PSYCHIATRIC RULES
        if any(g.startswith(('HTR', 'DRD')) for g in gene_symbols):
            return 'Psychiatric', 'Neurotransmitter Modulator'

        return 'Other', 'Experimental'

# Singleton Instance
drug_category_service = DrugCategoryService()
