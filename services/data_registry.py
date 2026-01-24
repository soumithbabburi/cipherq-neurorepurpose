"""
Unified Data Registry Service
Aggregates FDA drugs, chemical compounds, and provides consistent search interface
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from functools import lru_cache

logger = logging.getLogger(__name__)

class DataRegistry:
    """
    Unified data registry for drugs and compounds
    Provides single interface for all therapeutic category searches
    """
    
    THERAPEUTIC_CATEGORIES = [
        'cardiovascular',
        'diabetes', 
        'anti_inflammatory',
        'neurological',
        'psychiatric',
        'antibiotic',
        'antiviral',
        'cancer',
        'pain'
    ]
    
    CATEGORY_ALIASES = {
        'diabetic': 'diabetes',
        'heart': 'cardiovascular',
        'cardiac': 'cardiovascular',
        'inflammation': 'anti_inflammatory',
        'brain': 'neurological',
        'neuro': 'neurological',
        'mental': 'psychiatric',
        'psych': 'psychiatric',
        'depression': 'psychiatric',
        'anxiety': 'psychiatric',
        'infection': 'antibiotic',
        'bacterial': 'antibiotic',
        'virus': 'antiviral',
        'viral': 'antiviral',
        'oncology': 'cancer',
        'tumor': 'cancer',
        'analgesic': 'pain',
        'painkiller': 'pain',
    }
    
    COMPOUND_CATEGORIES = [
        'natural_products',
        'kinase_inhibitors',
        'gpcr_modulators',
        'protease_inhibitors',
        'epigenetic_modulators',
        'ion_channel_modulators',
        'metabolic_targets',
        'protein_protein_inhibitors',
        'degraders_protacs',
        'fragment_library'
    ]
    
    def __init__(self):
        self._drugs_40k = None
        self._compounds_10k = None
        self._drugs_by_category = {}
        self._compounds_by_category = {}
        self._loaded = False
        
    def _ensure_loaded(self):
        """Lazy load data on first access"""
        if self._loaded:
            return
        
        self._load_drugs()
        self._load_compounds()
        self._loaded = True
        
    def _load_drugs(self):
        """Load 40k drugs dataset"""
        try:
            drugs_path = Path(__file__).parent.parent / 'data' / 'drugs_40k.json'
            
            if drugs_path.exists():
                with open(drugs_path, 'r') as f:
                    self._drugs_40k = json.load(f)
                    
                for drug in self._drugs_40k:
                    category = drug.get('therapeutic_category', 'other')
                    if category not in self._drugs_by_category:
                        self._drugs_by_category[category] = []
                    self._drugs_by_category[category].append(drug)
                    
                logger.info(f"Loaded {len(self._drugs_40k)} drugs from 40k dataset")
            else:
                drugs_path_500 = Path(__file__).parent.parent / 'data' / 'drugs_500.json'
                if drugs_path_500.exists():
                    with open(drugs_path_500, 'r') as f:
                        self._drugs_40k = json.load(f)
                    logger.info(f"Loaded {len(self._drugs_40k)} drugs from 500 dataset (fallback)")
                else:
                    self._drugs_40k = []
                    logger.warning("No drugs dataset found")
                    
        except Exception as e:
            logger.error(f"Error loading drugs: {e}")
            self._drugs_40k = []
            
    def _load_compounds(self):
        """Load 10k compounds dataset"""
        try:
            compounds_path = Path(__file__).parent.parent / 'data' / 'compounds_10k.json'
            
            if compounds_path.exists():
                with open(compounds_path, 'r') as f:
                    self._compounds_10k = json.load(f)
                    
                for compound in self._compounds_10k:
                    category = compound.get('category', 'other')
                    if category not in self._compounds_by_category:
                        self._compounds_by_category[category] = []
                    self._compounds_by_category[category].append(compound)
                    
                logger.info(f"Loaded {len(self._compounds_10k)} compounds from 10k dataset")
            else:
                self._compounds_10k = []
                logger.warning("No compounds dataset found")
                
        except Exception as e:
            logger.error(f"Error loading compounds: {e}")
            self._compounds_10k = []
    
    def normalize_category(self, category: str) -> Optional[str]:
        """Normalize category name to standard form"""
        if not category:
            return None
            
        normalized = category.lower().strip().replace(' ', '_').replace('-', '_')
        
        if normalized in self.CATEGORY_ALIASES:
            return self.CATEGORY_ALIASES[normalized]
        
        if normalized in self.THERAPEUTIC_CATEGORIES:
            return normalized
            
        for cat in self.THERAPEUTIC_CATEGORIES:
            if normalized in cat or cat in normalized:
                return cat
                
        return normalized
    
    def search_drugs(
        self, 
        query: str = None,
        category: str = None,
        limit: int = 50,
        include_variants: bool = True
    ) -> List[Dict]:
        """
        Search drugs by query and/or category
        
        Args:
            query: Search term (name, class, target)
            category: Therapeutic category filter
            limit: Maximum results
            include_variants: Include drug variants (formulations)
            
        Returns:
            List of matching drugs
        """
        self._ensure_loaded()
        
        results = []
        
        normalized_category = self.normalize_category(category) if category else None
        
        if normalized_category and normalized_category in self._drugs_by_category:
            source_drugs = self._drugs_by_category[normalized_category]
        else:
            source_drugs = self._drugs_40k or []
        
        query_lower = query.lower() if query else None
        
        for drug in source_drugs:
            if query_lower:
                name = drug.get('name', '').lower()
                drug_class = drug.get('class', '').lower()
                target = drug.get('target', '').lower()
                
                if not (query_lower in name or 
                        query_lower in drug_class or 
                        query_lower in target):
                    continue
            
            if not include_variants:
                source = drug.get('source', '')
                if source == 'generated':
                    continue
            
            results.append(drug)
            
            if len(results) >= limit:
                break
        
        return results
    
    def search_compounds(
        self,
        query: str = None,
        category: str = None,
        limit: int = 50
    ) -> List[Dict]:
        """
        Search compounds by query and/or category
        
        Args:
            query: Search term (name, target, formula)
            category: Compound category filter
            limit: Maximum results
            
        Returns:
            List of matching compounds
        """
        self._ensure_loaded()
        
        results = []
        
        if category and category.lower() in self._compounds_by_category:
            source_compounds = self._compounds_by_category[category.lower()]
        else:
            source_compounds = self._compounds_10k or []
        
        query_lower = query.lower() if query else None
        
        for compound in source_compounds:
            if query_lower:
                name = compound.get('name', '').lower()
                target = compound.get('target', '').lower()
                formula = compound.get('formula', '').lower()
                
                if not (query_lower in name or 
                        query_lower in target or 
                        query_lower in formula):
                    continue
            
            results.append(compound)
            
            if len(results) >= limit:
                break
        
        return results
    
    def get_drugs_by_category(self, category: str, limit: int = 100) -> List[Dict]:
        """Get all drugs in a therapeutic category"""
        return self.search_drugs(category=category, limit=limit)
    
    def get_compounds_by_category(self, category: str, limit: int = 100) -> List[Dict]:
        """Get all compounds in a category"""
        return self.search_compounds(category=category, limit=limit)
    
    def get_drug_by_name(self, name: str) -> Optional[Dict]:
        """Get a specific drug by exact or partial name match"""
        results = self.search_drugs(query=name, limit=1)
        return results[0] if results else None
    
    def get_compound_by_name(self, name: str) -> Optional[Dict]:
        """Get a specific compound by name"""
        results = self.search_compounds(query=name, limit=1)
        return results[0] if results else None
    
    def get_drug_count(self) -> int:
        """Get total number of drugs"""
        self._ensure_loaded()
        return len(self._drugs_40k) if self._drugs_40k else 0
    
    def get_compound_count(self) -> int:
        """Get total number of compounds"""
        self._ensure_loaded()
        return len(self._compounds_10k) if self._compounds_10k else 0
    
    def get_drug_category_counts(self) -> Dict[str, int]:
        """Get drug counts by therapeutic category"""
        self._ensure_loaded()
        return {cat: len(drugs) for cat, drugs in self._drugs_by_category.items()}
    
    def get_compound_category_counts(self) -> Dict[str, int]:
        """Get compound counts by category"""
        self._ensure_loaded()
        return {cat: len(compounds) for cat, compounds in self._compounds_by_category.items()}
    
    def get_all_drug_categories(self) -> List[str]:
        """Get all available drug therapeutic categories"""
        return self.THERAPEUTIC_CATEGORIES.copy()
    
    def get_all_compound_categories(self) -> List[str]:
        """Get all available compound categories"""
        return self.COMPOUND_CATEGORIES.copy()
    
    def unified_search(
        self,
        query: str,
        search_type: str = 'all',
        category: str = None,
        limit: int = 50
    ) -> Dict[str, List[Dict]]:
        """
        Unified search across drugs and compounds
        
        Args:
            query: Search term
            search_type: 'drugs', 'compounds', or 'all'
            category: Category filter
            limit: Maximum results per type
            
        Returns:
            Dict with 'drugs' and 'compounds' lists
        """
        results = {'drugs': [], 'compounds': []}
        
        if search_type in ('drugs', 'all'):
            results['drugs'] = self.search_drugs(
                query=query, 
                category=category, 
                limit=limit
            )
        
        if search_type in ('compounds', 'all'):
            results['compounds'] = self.search_compounds(
                query=query,
                category=category,
                limit=limit
            )
        
        return results
    
    def get_statistics(self) -> Dict:
        """Get overall statistics about the data registry"""
        self._ensure_loaded()
        
        return {
            'total_drugs': self.get_drug_count(),
            'total_compounds': self.get_compound_count(),
            'drug_categories': len(self._drugs_by_category),
            'compound_categories': len(self._compounds_by_category),
            'drugs_by_category': self.get_drug_category_counts(),
            'compounds_by_category': self.get_compound_category_counts()
        }


_registry_instance = None

def get_data_registry() -> DataRegistry:
    """Get singleton data registry instance"""
    global _registry_instance
    if _registry_instance is None:
        _registry_instance = DataRegistry()
    return _registry_instance
