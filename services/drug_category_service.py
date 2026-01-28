import pandas as pd
import json
import logging
from typing import Dict, List, Tuple, Optional

logger = logging.getLogger(__name__)

class DrugCategoryService:
    def __init__(self):
        self.therapeutic_taxonomy = {
            'Cardiovascular': {'relevance_score': 0.72},
            'Metabolic': {'relevance_score': 0.68},
            'Anti-inflammatory': {'relevance_score': 0.65},
            'Neuroprotective': {'relevance_score': 0.78},
            'Psychiatric': {'relevance_score': 0.61},
            'Antibiotic_Repurposed': {'relevance_score': 0.58},
            'Other': {'relevance_score': 0.45}
        }
        self.gene_metadata = {}
        self._load_tsv_metadata()

    def _load_tsv_metadata(self):
        """Loads the functional categories from the TSV file"""
        try:
            # Note: keeping the name exactly as 'categories (1).tsv' as requested
            df = pd.read_csv('categories (1).tsv', sep='\t')
            # Map Gene Symbol -> List of functional categories (e.g. 'DPP4' -> ['PROTEASE'])
            self.gene_metadata = df.groupby('name')['name-2'].apply(set).to_dict()
            logger.info(f"Successfully mapped {len(self.gene_metadata)} genes from TSV")
        except Exception as e:
            logger.error(f"Error loading TSV: {e}")
            self.gene_metadata = {}

    def categorize_by_genes(self, gene_symbols: List[str]) -> Tuple[str, str]:
        """
        The core logic that solves the 'No Diabetic Drugs' issue.
        It uses the gene targets to determine the therapeutic category.
        """
        # Get functional classes from TSV for these genes
        functional_classes = set()
        for gene in gene_symbols:
            if gene in self.gene_metadata:
                functional_classes.update(self.gene_metadata[gene])

        # 1. METABOLIC LOGIC (Diabetes/Obesity)
        metabolic_markers = {'DPP4', 'PPARG', 'ABCC8', 'KCNJ11', 'SLC5A2', 'PRKAA1', 'PRKAA2', 'GLP1R'}
        if any(g in metabolic_markers for g in gene_symbols):
            return 'Metabolic', 'Diabetes/Insulin Regulator'

        # 2. CARDIOVASCULAR LOGIC
        if 'ION CHANNEL' in functional_classes or any(g.startswith(('ADR', 'ACE', 'KCN')) for g in gene_symbols):
            return 'Cardiovascular', 'Blood Pressure/Ion Channel'

        # 3. ANTI-INFLAMMATORY LOGIC
        if 'PROTEASE' in functional_classes or any(g in {'PTGS2', 'PTGS1', 'TNF', 'IL1B'} for g in gene_symbols):
            return 'Anti-inflammatory', 'NSAID/Inflammation'

        # 4. PSYCHIATRIC LOGIC
        if 'G PROTEIN COUPLED RECEPTOR' in functional_classes and any(g.startswith(('HTR', 'DRD')) for g in gene_symbols):
            return 'Psychiatric', 'Neurotransmitter Modulator'

        return 'Other', 'Experimental'

# Global instance
drug_category_service = DrugCategoryService()
