"""
Query Parser - Extract intent and filters from user queries
Handles category filtering, disease targets, and constraints
"""

import re
from typing import Dict, Optional, List
import logging

logger = logging.getLogger(__name__)

class QueryParser:
    """Parse user queries to extract category filters and intent"""
    
    def __init__(self):
        # Category keywords and their mappings
        self.category_keywords = {
            'cardiovascular': ['cardiovascular', 'cardiac', 'heart', 'statin', 'blood pressure', 
                              'ace inhibitor', 'beta blocker', 'anticoagulant'],
            'metabolic': ['metabolic', 'diabetes', 'insulin', 'glucose', 'metformin', 'ampk'],
            'anti-inflammatory': ['anti-inflammatory', 'inflammatory', 'inflammation', 
                                 'nsaid', 'cox inhibitor', 'corticosteroid'],
            'psychiatric': ['psychiatric', 'antidepressant', 'antipsychotic', 'ssri', 
                           'snri', 'anxiety', 'depression', 'mood'],
            'antibiotic_repurposed': ['antibiotic', 'antibacterial', 'penicillin', 
                                     'cephalosporin', 'antimicrobial'],
            'neuroprotective': ['neuroprotective', 'neuroprotection', 'antioxidant', 
                               'neurotrophic'],
            'antiviral_antiparasitic': ['antiviral', 'antiparasitic', 'antimalarial', 
                                       'antiviral'],
            'hormonal': ['hormonal', 'hormone', 'thyroid', 'estrogen', 'testosterone']
        }
        
        # Disease keywords
        self.disease_keywords = {
            'alzheimer': ["alzheimer", "alzheimer's", "ad", "dementia"],
            'parkinson': ["parkinson", "parkinson's", "pd"],
            'general': ['drug', 'repurposing', 'candidate']
        }
        
        # Constraint keywords
        self.constraint_keywords = {
            'fda_approved': ['fda approved', 'approved', 'licensed'],
            'in_trials': ['trial', 'clinical trial', 'in trials'],
            'safe': ['safe', 'safety', 'well-tolerated']
        }
    
    def parse_query(self, query: str) -> Dict:
        """
        Parse user query and extract structured information
        Returns: {
            'category_filter': str or None,
            'disease': str,
            'constraints': List[str],
            'original_query': str
        }
        """
        query_lower = query.lower()
        
        # Extract category filter
        category_filter = self._extract_category(query_lower)
        
        # Extract disease target
        disease = self._extract_disease(query_lower)
        
        # Extract constraints
        constraints = self._extract_constraints(query_lower)
        
        result = {
            'category_filter': category_filter,
            'disease': disease,
            'constraints': constraints,
            'original_query': query
        }
        
        logger.info(f"Parsed query: {result}")
        return result
    
    def _extract_category(self, query: str) -> Optional[str]:
        """Extract therapeutic category from query"""
        # Check for explicit "only" or "specific" keywords
        only_pattern = r'(\w+)\s+(drugs?\s+)?only'
        match = re.search(only_pattern, query)
        if match:
            category_term = match.group(1)
            for category, keywords in self.category_keywords.items():
                if category_term in keywords or category_term in category:
                    return self._normalize_category_name(category)
        
        # Check for category keywords in general
        for category, keywords in self.category_keywords.items():
            for keyword in keywords:
                if keyword in query:
                    # Check if this is truly a filter request
                    if any(term in query for term in ['only', 'just', 'specifically', 'focus on']):
                        return self._normalize_category_name(category)
        
        return None
    
    def _normalize_category_name(self, category: str) -> str:
        """Normalize category name to match DrugCategoryService"""
        mapping = {
            'cardiovascular': 'Cardiovascular',
            'metabolic': 'Metabolic',
            'anti-inflammatory': 'Anti-inflammatory',
            'psychiatric': 'Psychiatric',
            'antibiotic_repurposed': 'Antibiotic_Repurposed',
            'neuroprotective': 'Neuroprotective',
            'antiviral_antiparasitic': 'Antiviral_Antiparasitic',
            'hormonal': 'Hormonal'
        }
        return mapping.get(category, 'Other')
    
    def _extract_disease(self, query: str) -> str:
        """Extract disease target from query"""
        for disease, keywords in self.disease_keywords.items():
            for keyword in keywords:
                if keyword in query:
                    return disease
        return 'alzheimer'
    
    def _extract_constraints(self, query: str) -> List[str]:
        """Extract constraints from query"""
        constraints = []
        for constraint, keywords in self.constraint_keywords.items():
            for keyword in keywords:
                if keyword in query:
                    constraints.append(constraint)
                    break
        return constraints
    
    def format_filter_summary(self, parsed_query: Dict, drug_count: int) -> str:
        """Format a user-friendly summary of applied filters"""
        parts = []
        
        if parsed_query['category_filter']:
            parts.append(f"Category: {parsed_query['category_filter']}")
        else:
            parts.append("Category: All drugs")
        
        parts.append(f"Disease: {parsed_query['disease'].title()}")
        
        if parsed_query['constraints']:
            parts.append(f"Filters: {', '.join(parsed_query['constraints'])}")
        
        parts.append(f"Analyzing: {drug_count} drugs")
        
        return " | ".join(parts)


# Global instance
query_parser = QueryParser()
