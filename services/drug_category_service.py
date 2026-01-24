"""
Drug Category Service - Therapeutic Classification System
Organizes 500 drugs into therapeutic categories and subcategories for semantic reasoning
"""

import json
import logging
from typing import Dict, List, Set, Tuple, Optional
from pathlib import Path

logger = logging.getLogger(__name__)

class DrugCategoryService:
    """
    Comprehensive drug categorization system for repurposing analysis
    Maps drugs to therapeutic categories and Alzheimer's mechanisms
    """
    
    def __init__(self):
        # Major therapeutic categories with subcategories
        self.therapeutic_taxonomy = {
            'Cardiovascular': {
                'subcategories': ['Statins', 'ACE Inhibitors', 'Beta Blockers', 
                                 'Calcium Channel Blockers', 'ARBs', 'Anticoagulants',
                                 'Antiplatelet Agents', 'Diuretics'],
                'alzheimer_mechanisms': [
                    'Improves cerebral blood flow',
                    'Reduces vascular inflammation',
                    'Lowers cholesterol (reduces Aβ production)',
                    'Regulates blood pressure (reduces tau pathology)',
                    'Prevents microstrokes'
                ],
                'relevance_score': 0.72
            },
            'Metabolic': {
                'subcategories': ['Diabetes Drugs', 'AMPK Activators', 'PPARγ Agonists',
                                 'DPP4 Inhibitors', 'Metformin-like', 'Insulin Sensitizers'],
                'alzheimer_mechanisms': [
                    'Improves insulin signaling in brain',
                    'Reduces insulin resistance (linked to AD)',
                    'Activates AMPK (autophagy, energy metabolism)',
                    'Anti-inflammatory effects',
                    'Reduces AGE formation'
                ],
                'relevance_score': 0.68
            },
            'Anti-inflammatory': {
                'subcategories': ['NSAIDs', 'COX Inhibitors', 'Corticosteroids',
                                 'Immunomodulators', 'Anti-TNF Agents'],
                'alzheimer_mechanisms': [
                    'Reduces neuroinflammation',
                    'Lowers microglial activation',
                    'Decreases cytokine production',
                    'Protects blood-brain barrier',
                    'Reduces oxidative stress'
                ],
                'relevance_score': 0.65
            },
            'Neuroprotective': {
                'subcategories': ['Antioxidants', 'Neuropeptides', 'Neurotrophic Factors',
                                 'Calcium Modulators', 'Mitochondrial Protectants'],
                'alzheimer_mechanisms': [
                    'Direct neuroprotection',
                    'Reduces oxidative damage',
                    'Enhances mitochondrial function',
                    'Promotes neuronal survival',
                    'Supports synaptic plasticity'
                ],
                'relevance_score': 0.78
            },
            'Antibiotic_Repurposed': {
                'subcategories': ['Beta-lactams', 'Macrolides', 'Tetracyclines',
                                 'Fluoroquinolones', 'Antifungals'],
                'alzheimer_mechanisms': [
                    'Antimicrobial hypothesis (infections trigger AD)',
                    'Anti-inflammatory properties',
                    'Disrupts Aβ aggregation',
                    'Neuroprotective side effects',
                    'Microbiome modulation'
                ],
                'relevance_score': 0.58
            },
            'Psychiatric': {
                'subcategories': ['Antidepressants', 'Antipsychotics', 'Anxiolytics',
                                 'Mood Stabilizers', 'ADHD Medications'],
                'alzheimer_mechanisms': [
                    'Improves neuroplasticity',
                    'Enhances neurotransmitter balance',
                    'Reduces stress-induced damage',
                    'Protects against depression-related cognitive decline',
                    'Modulates serotonin/dopamine pathways'
                ],
                'relevance_score': 0.61
            },
            'Hormonal': {
                'subcategories': ['Hormone Replacement', 'Thyroid Drugs', 
                                 'Corticosteroids', 'Sex Hormones', 'Growth Factors'],
                'alzheimer_mechanisms': [
                    'Hormonal regulation of cognition',
                    'Neuroprotective hormone effects',
                    'Reduces inflammation',
                    'Supports synaptic function',
                    'Modulates Aβ clearance'
                ],
                'relevance_score': 0.63
            },
            'Antiviral_Antiparasitic': {
                'subcategories': ['Antivirals', 'Antiparasitics', 'Antimalarials',
                                 'Antiretrovirals', 'Anthelmintics'],
                'alzheimer_mechanisms': [
                    'Infectious hypothesis (HSV-1, others linked to AD)',
                    'Immunomodulation',
                    'Anti-inflammatory',
                    'Autophagy activation',
                    'Neuroprotection'
                ],
                'relevance_score': 0.56
            },
            'Other': {
                'subcategories': ['Supplements', 'Vitamins', 'Herbal', 'Experimental'],
                'alzheimer_mechanisms': [
                    'Various mechanisms',
                    'Nutritional support',
                    'Antioxidant effects'
                ],
                'relevance_score': 0.45
            }
        }
        
        # Drug class to category mapping
        self.class_to_category_map = {
            # Cardiovascular
            'Statin': ('Cardiovascular', 'Statins'),
            'HMG-CoA Reductase Inhibitor': ('Cardiovascular', 'Statins'),
            'ACE Inhibitor': ('Cardiovascular', 'ACE Inhibitors'),
            'Beta Blocker': ('Cardiovascular', 'Beta Blockers'),
            'Beta-Blocker': ('Cardiovascular', 'Beta Blockers'),
            'Calcium Channel Blocker': ('Cardiovascular', 'Calcium Channel Blockers'),
            'ARB': ('Cardiovascular', 'ARBs'),
            'Angiotensin II Receptor Blocker': ('Cardiovascular', 'ARBs'),
            'Anticoagulant': ('Cardiovascular', 'Anticoagulants'),
            'Anticoagulants': ('Cardiovascular', 'Anticoagulants'),
            'Antiplatelet': ('Cardiovascular', 'Antiplatelet Agents'),
            'Diuretic': ('Cardiovascular', 'Diuretics'),
            'Antihypertensives': ('Cardiovascular', 'ACE Inhibitors'),
            'Antiarrhythmic': ('Cardiovascular', 'Antiarrhythmics'),
            'Alpha Blocker': ('Cardiovascular', 'Alpha Blockers'),
            'Bile Acid Sequestrant': ('Cardiovascular', 'Statins'),
            
            # Metabolic
            'Diabetes Drug': ('Metabolic', 'Diabetes Drugs'),
            'Antidiabetic': ('Metabolic', 'Diabetes Drugs'),
            'Biguanide': ('Metabolic', 'Metformin-like'),
            'Thiazolidinedione': ('Metabolic', 'PPARγ Agonists'),
            'DPP4 Inhibitor': ('Metabolic', 'DPP4 Inhibitors'),
            'DPP-4 Inhibitor': ('Metabolic', 'DPP4 Inhibitors'),
            'SGLT2 Inhibitor': ('Metabolic', 'Diabetes Drugs'),
            'GLP-1 Agonist': ('Metabolic', 'Diabetes Drugs'),
            'Sulfonylurea': ('Metabolic', 'Diabetes Drugs'),
            'Rapid-Acting Insulin': ('Metabolic', 'Insulin Sensitizers'),
            'AMPK Activator': ('Metabolic', 'AMPK Activators'),
            'Insulin Sensitizer': ('Metabolic', 'Insulin Sensitizers'),
            
            # Anti-inflammatory
            'NSAID': ('Anti-inflammatory', 'NSAIDs'),
            'COX Inhibitor': ('Anti-inflammatory', 'COX Inhibitors'),
            'COX-2 Inhibitor': ('Anti-inflammatory', 'COX Inhibitors'),
            'Corticosteroid': ('Anti-inflammatory', 'Corticosteroids'),
            'Immunomodulator': ('Anti-inflammatory', 'Immunomodulators'),
            'Anti-TNF': ('Anti-inflammatory', 'Anti-TNF Agents'),
            'TNF Inhibitor': ('Anti-inflammatory', 'Anti-TNF Agents'),
            'Immunosuppressants': ('Anti-inflammatory', 'Immunomodulators'),
            'JAK Inhibitor': ('Anti-inflammatory', 'Immunomodulators'),
            'Bisphosphonate': ('Anti-inflammatory', 'Immunomodulators'),
            
            # Antibiotics
            'Antibiotic': ('Antibiotic_Repurposed', 'Beta-lactams'),
            'Beta-lactam': ('Antibiotic_Repurposed', 'Beta-lactams'),
            'Penicillin': ('Antibiotic_Repurposed', 'Beta-lactams'),
            'Cephalosporin': ('Antibiotic_Repurposed', 'Beta-lactams'),
            'Macrolide': ('Antibiotic_Repurposed', 'Macrolides'),
            'Tetracycline': ('Antibiotic_Repurposed', 'Tetracyclines'),
            'Fluoroquinolone': ('Antibiotic_Repurposed', 'Fluoroquinolones'),
            'Antifungal': ('Antibiotic_Repurposed', 'Antifungals'),
            'Antifungals': ('Antibiotic_Repurposed', 'Antifungals'),
            
            # Psychiatric
            'Antidepressant': ('Psychiatric', 'Antidepressants'),
            'SSRI': ('Psychiatric', 'Antidepressants'),
            'SNRI': ('Psychiatric', 'Antidepressants'),
            'TCA': ('Psychiatric', 'Antidepressants'),
            'MAOI': ('Psychiatric', 'Antidepressants'),
            'Antipsychotic': ('Psychiatric', 'Antipsychotics'),
            'Atypical Antipsychotic': ('Psychiatric', 'Antipsychotics'),
            'Anxiolytic': ('Psychiatric', 'Anxiolytics'),
            'Mood Stabilizer': ('Psychiatric', 'Mood Stabilizers'),
            'Anticonvulsant': ('Psychiatric', 'Mood Stabilizers'),
            
            # Neuroprotective
            'Antioxidant': ('Neuroprotective', 'Antioxidants'),
            'Neuroprotective': ('Neuroprotective', 'Neuropeptides'),
            'Neurotrophic': ('Neuroprotective', 'Neurotrophic Factors'),
            'Cholinesterase Inhibitor': ('Neuroprotective', 'Neuropeptides'),
            'Dopamine Agonist': ('Neuroprotective', 'Neuropeptides'),
            'Retinoid': ('Neuroprotective', 'Antioxidants'),
            
            # Antivirals
            'Antiviral': ('Antiviral_Antiparasitic', 'Antivirals'),
            'Antivirals': ('Antiviral_Antiparasitic', 'Antivirals'),
            'Antiparasitic': ('Antiviral_Antiparasitic', 'Antiparasitics'),
            'Antimalarial': ('Antiviral_Antiparasitic', 'Antimalarials'),
            
            # Hormonal
            'Hormone': ('Hormonal', 'Hormone Replacement'),
            'Hormones': ('Hormonal', 'Hormone Replacement'),
            'Thyroid Drug': ('Hormonal', 'Thyroid Drugs'),
            'Estrogen': ('Hormonal', 'Sex Hormones'),
            'Testosterone': ('Hormonal', 'Sex Hormones'),
            
            # Other Respiratory/Cancer/etc categories that are less relevant to AD
            # But still categorize them for completeness
            'Bronchodilators': ('Other', 'Respiratory'),
            'LABA': ('Other', 'Respiratory'),
            'LAMA': ('Other', 'Respiratory'),
            'Antihistamines': ('Other', 'Allergy'),
            '5-HT3 Antagonist': ('Other', 'Allergy'),
            'Muscle Relaxants': ('Other', 'Musculoskeletal'),
            'Opioid': ('Other', 'Pain'),
            'Aminoglycoside': ('Antibiotic_Repurposed', 'Beta-lactams'),
            
            # Cancer drugs - less relevant for AD but categorize for completeness
            'EGFR Inhibitor': ('Other', 'Cancer'),
            'Chemotherapy': ('Other', 'Cancer'),
            'TKI': ('Other', 'Cancer'),
            'ALK Inhibitor': ('Other', 'Cancer'),
            'Multi-Kinase Inhibitor': ('Other', 'Cancer'),
            'HER2 Inhibitor': ('Other', 'Cancer'),
            'PARP Inhibitor': ('Other', 'Cancer'),
            'Aromatase Inhibitor': ('Other', 'Cancer'),
            'Antiandrogen': ('Other', 'Cancer'),
            'Antimetabolite': ('Other', 'Cancer'),
            'Topoisomerase Inhibitor': ('Other', 'Cancer'),
        }
    
    def categorize_drug(self, drug: Dict) -> Tuple[str, str]:
        """
        Categorize a single drug into therapeutic category and subcategory
        Returns: (category, subcategory)
        """
        drug_class = drug.get('class', 'Unknown')
        
        # Try exact match first
        if drug_class in self.class_to_category_map:
            return self.class_to_category_map[drug_class]
        
        # Try partial match
        for class_name, (category, subcategory) in self.class_to_category_map.items():
            if class_name.lower() in drug_class.lower() or drug_class.lower() in class_name.lower():
                return (category, subcategory)
        
        # Default to Other
        return ('Other', 'Experimental')
    
    def categorize_all_drugs(self, drugs: List[Dict]) -> Dict[str, Dict]:
        """
        Categorize all drugs and organize by therapeutic category
        Returns nested dict: {category: {subcategory: [drugs]}}
        """
        categorized = {}
        
        for drug in drugs:
            category, subcategory = self.categorize_drug(drug)
            
            if category not in categorized:
                categorized[category] = {}
            
            if subcategory not in categorized[category]:
                categorized[category][subcategory] = []
            
            # Add category info to drug
            drug_with_category = drug.copy()
            drug_with_category['therapeutic_category'] = category
            drug_with_category['subcategory'] = subcategory
            
            categorized[category][subcategory].append(drug_with_category)
        
        logger.info(f"Categorized {len(drugs)} drugs into {len(categorized)} therapeutic categories")
        return categorized
    
    def get_category_info(self, category: str) -> Optional[Dict]:
        """Get information about a therapeutic category"""
        return self.therapeutic_taxonomy.get(category)
    
    def filter_by_category(self, drugs: List[Dict], category_filter: Optional[str] = None) -> List[Dict]:
        """
        Filter drugs by therapeutic category
        If category_filter is None, return all drugs
        """
        if not category_filter:
            return drugs
        
        categorized = self.categorize_all_drugs(drugs)
        
        if category_filter not in categorized:
            logger.warning(f"Category '{category_filter}' not found. Available: {list(categorized.keys())}")
            return []
        
        # Flatten all drugs in this category
        filtered_drugs = []
        for subcategory, drug_list in categorized[category_filter].items():
            filtered_drugs.extend(drug_list)
        
        logger.info(f"Filtered to {len(filtered_drugs)} drugs in category '{category_filter}'")
        return filtered_drugs
    
    def get_category_statistics(self, drugs: List[Dict]) -> Dict:
        """Get statistics about drug distribution across categories"""
        categorized = self.categorize_all_drugs(drugs)
        
        stats = {
            'total_drugs': len(drugs),
            'categories': {}
        }
        
        for category, subcategories in categorized.items():
            total_in_category = sum(len(drugs) for drugs in subcategories.values())
            stats['categories'][category] = {
                'total_drugs': total_in_category,
                'subcategories': {subcat: len(drugs) for subcat, drugs in subcategories.items()},
                'relevance_score': self.therapeutic_taxonomy.get(category, {}).get('relevance_score', 0.5),
                'mechanisms': self.therapeutic_taxonomy.get(category, {}).get('alzheimer_mechanisms', [])
            }
        
        return stats
    
    def get_top_drugs_per_category(self, scored_drugs: List[Dict], 
                                   top_n: int = 3, 
                                   category_filter: Optional[str] = None) -> Dict:
        """
        Get top N drugs from each therapeutic category (or specific category if filtered)
        Input: scored_drugs = [{'name': 'Drug', 'score': 0.85, 'class': '...', ...}, ...]
        Returns: {category: {subcategory: [top_drugs]}}
        """
        # First categorize all drugs
        categorized = self.categorize_all_drugs(scored_drugs)
        
        # Filter to specific category if requested
        if category_filter and category_filter in categorized:
            categorized = {category_filter: categorized[category_filter]}
        
        # Get top N from each subcategory
        top_drugs_by_category = {}
        
        for category, subcategories in categorized.items():
            top_drugs_by_category[category] = {}
            
            for subcategory, drugs in subcategories.items():
                # Sort by score descending
                sorted_drugs = sorted(drugs, key=lambda x: x.get('score', 0), reverse=True)
                top_drugs_by_category[category][subcategory] = sorted_drugs[:top_n]
        
        return top_drugs_by_category


# Global instance
drug_category_service = DrugCategoryService()
