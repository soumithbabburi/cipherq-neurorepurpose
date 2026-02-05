"""
Intelligent Query Processor - ALL-IN-ONE
Understands: categories, proteins, pathways, diseases - using YOUR data only!
"""

import json
import logging
from typing import List, Dict

logger = logging.getLogger(__name__)


class IntelligentQueryProcessor:
    """Semantic query understanding using all platform data"""
    
    def __init__(self):
        # Load all data
        with open('drugs.json', 'r') as f:
            self.drugs = json.load(f)
        
        with open('drug_interactions.json', 'r') as f:
            self.interactions = json.load(f)
        
        with open('protein_pathways.json', 'r') as f:
            self.protein_pathways = json.load(f)
        
        with open('pathways.json', 'r') as f:
            self.pathways = json.load(f)
        
        with open('drug_therapeutic_categories.json', 'r') as f:
            self.categories = json.load(f)
        
        with open('genes.json', 'r') as f:
            self.genes = json.load(f)
        
        logger.info("✅ Intelligent query processor initialized")
    
    def process_query(self, query: str, target_disease: str = None) -> Dict:
        """
        Process natural language query and return relevant drugs
        
        Understands:
        - "dopaminergic drugs"
        - "drugs targeting PPARG"
        - "drugs in insulin signaling pathway"
        - "drugs for diabetes connected to Alzheimer's"
        """
        
        query_lower = query.lower()
        results = {'drugs': [], 'interpretation': '', 'method': ''}
        
        # PATTERN 1: Category queries (dopaminergic, diabetic, cardiovascular)
        if any(kw in query_lower for kw in ['dopaminergic', 'cholinergic', 'diabetic', 'cardiovascular', 'pain', 'psychiatric']):
            results = self._process_category_query(query_lower)
        
        # PATTERN 2: Protein target queries (targeting PPARG, binds to DRD2)
        elif any(kw in query_lower for kw in ['targeting', 'target', 'binds', 'inhibit', 'agonist', 'antagonist']):
            results = self._process_protein_query(query_lower)
        
        # PATTERN 3: Pathway queries (insulin signaling, dopamine pathway)
        elif any(kw in query_lower for kw in ['pathway', 'signaling', 'metabolism']):
            results = self._process_pathway_query(query_lower)
        
        # PATTERN 4: Disease connection queries
        elif any(kw in query_lower for kw in ['for', 'treat', 'disease']):
            results = self._process_disease_query(query_lower, target_disease)
        
        # PATTERN 5: General drug search
        else:
            results = self._process_general_query(query_lower)
        
        return results
    
    def _process_category_query(self, query: str) -> Dict:
        """Handle category queries: 'dopaminergic drugs'"""
        
        category_map = {
            'dopaminergic': 'Parkinsons',
            'cholinergic': 'Alzheimers',
            'diabetic': 'Diabetic',
            'metabolic': 'Diabetic',
            'cardiovascular': 'Cardiovascular',
            'cardiac': 'Cardiovascular',
            'pain': 'Pain',
            'psychiatric': 'Psychiatric'
        }
        
        matched_category = None
        for keyword, category in category_map.items():
            if keyword in query:
                matched_category = category
                break
        
        if matched_category:
            drugs = self._get_drugs_by_category(matched_category)
            return {
                'drugs': drugs,
                'interpretation': f"Found {len(drugs)} drugs in {matched_category} category",
                'method': 'category_match'
            }
        
        return {'drugs': [], 'interpretation': 'No category matched', 'method': 'none'}
    
    def _process_protein_query(self, query: str) -> Dict:
        """Handle protein queries: 'drugs targeting PPARG'"""
        
        # Extract protein name from query
        # Look for uppercase words (gene symbols)
        words = query.split()
        protein_candidates = [w.upper() for w in words if len(w) >= 3 and w.isupper()]
        
        # Also check for common protein names
        common_proteins = {
            'pparg': 'PPARG', 'ppar': 'PPARG',
            'drd2': 'DRD2', 'dopamine receptor': 'DRD2',
            'ache': 'ACHE', 'acetylcholinesterase': 'ACHE',
            'dpp4': 'DPP4', 'ampk': 'PRKAA1',
            'cox2': 'PTGS2', 'cyclooxygenase': 'PTGS2'
        }
        
        for keyword, gene in common_proteins.items():
            if keyword in query:
                protein_candidates.append(gene)
        
        # Find drugs targeting these proteins
        matched_drugs = []
        
        for protein in protein_candidates:
            for drug_name, targets in self.interactions.items():
                for target in targets:
                    if target.get('gene_symbol', '').upper() == protein:
                        if drug_name not in [d['name'] for d in matched_drugs]:
                            matched_drugs.append({
                                'name': drug_name,
                                'target': protein,
                                'confidence': target.get('confidence_score', 0)
                            })
        
        if matched_drugs:
            return {
                'drugs': matched_drugs,
                'interpretation': f"Found {len(matched_drugs)} drugs targeting {', '.join(protein_candidates[:3])}",
                'method': 'protein_target'
            }
        
        return {'drugs': [], 'interpretation': 'No protein matches found', 'method': 'none'}
    
    def _process_pathway_query(self, query: str) -> Dict:
        """Handle pathway queries: 'drugs in insulin signaling pathway'"""
        
        # Extract pathway keywords
        pathway_keywords = []
        if 'insulin' in query:
            pathway_keywords.append('insulin')
        if 'glucose' in query or 'metabol' in query:
            pathway_keywords.append('glucose')
        if 'dopamin' in query:
            pathway_keywords.append('dopamin')
        if 'ampk' in query:
            pathway_keywords.append('ampk')
        
        # Find pathways matching keywords
        matched_pathway_ids = []
        for pw_id, pw_data in self.pathways.items():
            pw_name = pw_data.get('name', '').lower()
            if any(kw in pw_name for kw in pathway_keywords):
                matched_pathway_ids.append(pw_id)
        
        # Find proteins in these pathways
        proteins_in_pathways = []
        for gene, pathway_list in self.protein_pathways.items():
            if any(pw_id in pathway_list for pw_id in matched_pathway_ids):
                proteins_in_pathways.append(gene)
        
        # Find drugs targeting these proteins
        matched_drugs = []
        for drug_name, targets in self.interactions.items():
            for target in targets:
                if target.get('gene_symbol', '').upper() in proteins_in_pathways:
                    if drug_name not in [d['name'] for d in matched_drugs]:
                        matched_drugs.append({
                            'name': drug_name,
                            'pathway_match': ', '.join(pathway_keywords),
                            'target': target.get('gene_symbol')
                        })
                        break
        
        if matched_drugs:
            return {
                'drugs': matched_drugs,
                'interpretation': f"Found {len(matched_drugs)} drugs in {', '.join(pathway_keywords)} pathways",
                'method': 'pathway_match'
            }
        
        return {'drugs': [], 'interpretation': 'No pathway matches', 'method': 'none'}
    
    def _process_disease_query(self, query: str, target_disease: str) -> Dict:
        """Handle disease queries: 'drugs for Alzheimer's'"""
        
        # Use target_disease if provided, otherwise extract from query
        disease = target_disease or query
        
        # Get category most relevant to disease
        if 'alzheimer' in disease.lower() or 'dementia' in disease.lower():
            category = 'Diabetic'  # Metabolic drugs often relevant
        elif 'parkinson' in disease.lower():
            category = 'Parkinsons'
        elif 'diabetes' in disease.lower():
            category = 'Diabetic'
        else:
            category = 'General'
        
        drugs = self._get_drugs_by_category(category)
        
        return {
            'drugs': drugs,
            'interpretation': f"Showing {category} drugs potentially relevant to {disease}",
            'method': 'disease_category'
        }
    
    def _process_general_query(self, query: str) -> Dict:
        """Handle general searches"""
        
        # Search drug names
        matches = []
        for drug_name in self.drugs.keys():
            if query in drug_name:
                matches.append({'name': drug_name})
        
        return {
            'drugs': matches,
            'interpretation': f"Drug name search results",
            'method': 'name_search'
        }
    
    def _get_drugs_by_category(self, category: str) -> List[Dict]:
        """Get all drugs in a category"""
        
        drugs_in_category = []
        
        for drug_name, drug_categories in self.categories.items():
            if category in drug_categories:
                # Get targets
                targets = []
                if drug_name in self.interactions:
                    targets = [t.get('gene_symbol') for t in self.interactions[drug_name][:5]]
                
                drugs_in_category.append({
                    'name': drug_name,
                    'category': category,
                    'targets': targets
                })
        
        return drugs_in_category


def process_intelligent_query(query: str, target_disease: str = None) -> Dict:
    """
    Main entry point for intelligent query processing
    Call this from your app!
    """
    processor = IntelligentQueryProcessor()
    return processor.process_query(query, target_disease)


__all__ = ['IntelligentQueryProcessor', 'process_intelligent_query']
