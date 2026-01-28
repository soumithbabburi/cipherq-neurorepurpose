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
        try:
            # Robust path finding: assumes this file is in app/services/
            # and TSV is in app/
            current_dir = os.path.dirname(os.path.abspath(__file__))
            tsv_path = os.path.join(current_dir, '..', 'categories (1).tsv')
            tsv_path = os.path.normpath(tsv_path)
            
            if os.path.exists(tsv_path):
                df = pd.read_csv(tsv_path, sep='\t')
                self.gene_metadata = df.groupby('name')['name-2'].apply(set).to_dict()
                logger.info(f"✅ Loaded TSV from {tsv_path}")
            else:
                logger.warning(f"⚠️ TSV not found at {tsv_path}")
        except Exception as e:
            logger.error(f"❌ Error loading TSV: {e}")

    def categorize_by_genes(self, gene_symbols: List[str]) -> Tuple[str, str]:
        functions = set()
        for gene in gene_symbols:
            if gene in self.gene_metadata:
                functions.update(self.gene_metadata[gene])

        # --- ENHANCED DIABETES RULES ---
        # Added: INSR (Insulin Receptor), IGF1R, AKT, MTOR
        metabolic_markers = {
            'DPP4', 'PPARG', 'ABCC8', 'KCNJ11', 'SLC5A2', 'PRKAA1', 'PRKAA2', 
            'GLP1R', 'INSR', 'IGF1R', 'GCGR', 'SGLT2'
        }
        
        # Check for direct gene hits
        if any(g in metabolic_markers for g in gene_symbols):
            return 'Metabolic', 'Diabetes/Insulin Regulator'

        # Check for functional hits from TSV
        if 'INSULIN SIGNALING' in functions or 'GLUCOSE METABOLISM' in functions:
            return 'Metabolic', 'Diabetes/Insulin Regulator'

        # Cardiovascular
        if 'ION CHANNEL' in functions or any(g.startswith(('ADR', 'ACE')) for g in gene_symbols):
            return 'Cardiovascular', 'Vascular Support'

        # Anti-inflammatory
        if 'PROTEASE' in functions or any(g in {'PTGS2', 'TNF', 'IL1B'} for g in gene_symbols):
            return 'Anti-inflammatory', 'Neuroinflammation Control'

        return 'Other', 'Experimental'

drug_category_service = DrugCategoryService()
