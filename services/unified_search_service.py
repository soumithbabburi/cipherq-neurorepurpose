"""
Unified Search Service
Clean interface for searching both FDA drugs and chemical compounds
"""

import logging
from typing import List, Dict, Optional, Tuple
from enum import Enum

logger = logging.getLogger(__name__)


class SearchMode(Enum):
    """Search mode selection"""
    DRUGS = "drugs"  # FDA-approved drugs for repurposing
    COMPOUNDS = "compounds"  # Chemical compounds for discovery
    BOTH = "both"  # Search both databases


class UnifiedSearchService:
    """
    Unified service for searching drugs and compounds
    Provides clean interface for UI integration
    """
    
    # Therapeutic categories with display names
    DRUG_CATEGORIES = {
        'cardiovascular': 'Cardiovascular',
        'diabetes': 'Diabetes',
        'anti_inflammatory': 'Anti-inflammatory',
        'neurological': 'Neurological',
        'psychiatric': 'Psychiatric',
        'antibiotic': 'Antibiotic',
        'antiviral': 'Antiviral',
        'cancer': 'Cancer/Oncology',
        'pain': 'Pain Management',
        'gastrointestinal': 'Gastrointestinal',
        'respiratory': 'Respiratory',
        'immunology': 'Immunology',
        'endocrine': 'Endocrine',
        'other': 'Other'
    }
    
    def __init__(self):
        self._drug_service = None
        self._compound_service = None
    
    @property
    def drug_service(self):
        """Lazy load drug service"""
        if self._drug_service is None:
            from services.drug_database_service import get_drug_database_service
            self._drug_service = get_drug_database_service()
        return self._drug_service
    
    @property
    def compound_service(self):
        """Lazy load compound service"""
        if self._compound_service is None:
            from services.compound_database_service import get_compound_database_service
            self._compound_service = get_compound_database_service()
        return self._compound_service
    
    def search(
        self,
        query: str = None,
        mode: SearchMode = SearchMode.DRUGS,
        category: str = None,
        min_druglikeness: float = None,
        limit: int = 50
    ) -> Tuple[List[Dict], str]:
        """
        Unified search across drugs and/or compounds
        
        Args:
            query: Search term
            mode: Search mode (drugs, compounds, or both)
            category: Therapeutic category filter (drugs only)
            min_druglikeness: Minimum drug-likeness score (compounds only)
            limit: Maximum results
            
        Returns:
            Tuple of (results list, source description)
        """
        results = []
        source = ""
        
        if mode in (SearchMode.DRUGS, SearchMode.BOTH):
            drugs = self.drug_service.search_drugs(
                query=query,
                category=category,
                limit=limit
            )
            for drug in drugs:
                drug['_source'] = 'fda_drug'
            results.extend(drugs)
            source = "FDA-Approved Drugs"
        
        if mode in (SearchMode.COMPOUNDS, SearchMode.BOTH):
            remaining_limit = limit - len(results)
            if remaining_limit > 0:
                compounds = self.compound_service.search_compounds(
                    query=query,
                    min_druglikeness=min_druglikeness,
                    limit=remaining_limit
                )
                for compound in compounds:
                    compound['_source'] = 'chemical_compound'
                results.extend(compounds)
                
                if mode == SearchMode.COMPOUNDS:
                    source = "Chemical Compounds"
                else:
                    source = "FDA Drugs + Chemical Compounds"
        
        return results, source
    
    def get_drugs_by_category(
        self, 
        category: str, 
        limit: int = 100
    ) -> List[Dict]:
        """Get all drugs in a therapeutic category"""
        return self.drug_service.get_drugs_by_category(category, limit)
    
    def get_compound_by_smiles(self, smiles: str) -> Optional[Dict]:
        """Get or create compound from SMILES"""
        compound = self.compound_service.fetch_compound_by_smiles(smiles)
        if compound:
            self.compound_service.save_compound(compound)
        return compound
    
    def get_drug_by_name(self, name: str) -> Optional[Dict]:
        """Get drug by name"""
        return self.drug_service.get_drug_by_name(name)
    
    def get_database_stats(self) -> Dict:
        """Get statistics for both databases"""
        return {
            'fda_drugs': {
                'total': self.drug_service.get_drug_count(),
                'categories': self.drug_service.get_category_counts()
            },
            'compounds': {
                'total': self.compound_service.get_compound_count()
            }
        }
    
    def get_categories(self) -> List[Tuple[str, str]]:
        """Get list of therapeutic categories"""
        return list(self.DRUG_CATEGORIES.items())
    
    def quick_search(
        self, 
        query: str, 
        mode: str = "drugs"
    ) -> List[Dict]:
        """
        Quick search for autocomplete/instant results
        
        Args:
            query: Search term
            mode: "drugs", "compounds", or "both"
            
        Returns:
            List of matching items (max 10)
        """
        search_mode = SearchMode(mode) if mode in ['drugs', 'compounds', 'both'] else SearchMode.DRUGS
        results, _ = self.search(query=query, mode=search_mode, limit=10)
        return results


# Singleton instance
_unified_service = None

def get_unified_search_service() -> UnifiedSearchService:
    """Get singleton unified search service"""
    global _unified_service
    if _unified_service is None:
        _unified_service = UnifiedSearchService()
    return _unified_service
