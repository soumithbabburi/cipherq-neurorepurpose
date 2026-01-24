"""
Alzheimer's Drug Exclusion Filter - Remove existing AD drugs and their data
"""
import json
import logging
from typing import List, Dict, Any, Set
from pathlib import Path
import pandas as pd
from core.config_loader import config

logger = logging.getLogger(__name__)

class AlzheimerExclusionFilter:
    """Filter to exclude existing Alzheimer's drugs and their associated data"""
    
    def __init__(self):
        self.ad_drugs_data = config.load_alzheimer_drugs()
        self.excluded_drugs = self._build_exclusion_set()
        self.excluded_mesh_terms = set(self.ad_drugs_data.get("mesh_terms_to_exclude", []))
        
    def _build_exclusion_set(self) -> Set[str]:
        """Build set of drug names and synonyms to exclude"""
        exclusion_set = set()
        
        for drug_info in self.ad_drugs_data.get("existing_ad_drugs", []):
            if drug_info.get("exclude", True):
                # Add primary name
                exclusion_set.add(drug_info["name"].lower())
                
                # Add all synonyms
                for synonym in drug_info.get("synonyms", []):
                    exclusion_set.add(synonym.lower())
                    
                # Add DrugBank ID if present
                if "drugbank_id" in drug_info:
                    exclusion_set.add(drug_info["drugbank_id"].lower())
                    
        logger.info(f"Built Alzheimer's exclusion set with {len(exclusion_set)} terms")
        return exclusion_set
    
    def is_alzheimer_drug(self, drug_name: str) -> bool:
        """Check if a drug is already used for Alzheimer's disease"""
        if not drug_name:
            return False
            
        # Normalize drug name for comparison
        normalized_name = drug_name.lower().strip()
        
        # Check direct matches and common variations
        variations = [
            normalized_name,
            normalized_name.replace(" ", ""),
            normalized_name.replace("-", ""),
            normalized_name.replace("_", "")
        ]
        
        for variation in variations:
            if variation in self.excluded_drugs:
                return True
                
        return False
    
    def filter_drug_candidates(self, candidates_df: pd.DataFrame) -> pd.DataFrame:
        """Remove Alzheimer's drugs from candidate list"""
        if candidates_df.empty:
            return candidates_df
            
        # Determine the drug name column
        drug_col = None
        possible_cols = ["Drug", "drug", "name", "drug_name", "compound_name"]
        for col in possible_cols:
            if col in candidates_df.columns:
                drug_col = col
                break
                
        if not drug_col:
            logger.warning("No drug name column found in candidates DataFrame")
            return candidates_df
        
        # Filter out Alzheimer's drugs
        initial_count = len(candidates_df)
        filtered_df = candidates_df[~candidates_df[drug_col].apply(self.is_alzheimer_drug)].copy()
        filtered_count = len(filtered_df)
        
        excluded_count = initial_count - filtered_count
        if excluded_count > 0:
            logger.info(f"Excluded {excluded_count} existing Alzheimer's drugs from {initial_count} candidates")
            
        return filtered_df.reset_index(drop=True)
    
    def filter_publications(self, publications: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove publications related to existing Alzheimer's drugs"""
        if not publications:
            return publications
            
        filtered_publications = []
        excluded_count = 0
        
        for pub in publications:
            title = pub.get("title", "").lower()
            abstract = pub.get("abstract", "").lower()
            mesh_terms = pub.get("mesh_terms", [])
            
            # Check if publication is about excluded drugs
            exclude_pub = False
            
            # Check title and abstract for drug names
            for excluded_drug in self.excluded_drugs:
                if excluded_drug in title or excluded_drug in abstract:
                    exclude_pub = True
                    break
            
            # Check MeSH terms
            if not exclude_pub:
                for mesh_term in mesh_terms:
                    if mesh_term in self.excluded_mesh_terms:
                        exclude_pub = True
                        break
            
            if not exclude_pub:
                filtered_publications.append(pub)
            else:
                excluded_count += 1
                
        if excluded_count > 0:
            logger.info(f"Excluded {excluded_count} publications related to existing AD drugs")
            
        return filtered_publications
    
    def filter_clinical_trials(self, trials: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove clinical trials for existing Alzheimer's drugs"""
        if not trials:
            return trials
            
        filtered_trials = []
        excluded_count = 0
        
        for trial in trials:
            title = trial.get("title", "").lower()
            intervention = trial.get("intervention", "").lower()
            condition = trial.get("condition", "").lower()
            
            # Check if trial involves excluded drugs
            exclude_trial = False
            
            # Check for drug names in title, intervention, or condition
            for excluded_drug in self.excluded_drugs:
                if (excluded_drug in title or 
                    excluded_drug in intervention or 
                    excluded_drug in condition):
                    exclude_trial = True
                    break
            
            if not exclude_trial:
                filtered_trials.append(trial)
            else:
                excluded_count += 1
                
        if excluded_count > 0:
            logger.info(f"Excluded {excluded_count} clinical trials for existing AD drugs")
            
        return filtered_trials
    
    def get_excluded_drugs_info(self) -> Dict[str, Any]:
        """Get information about excluded Alzheimer's drugs"""
        return {
            "excluded_drugs": list(self.excluded_drugs),
            "excluded_count": len(self.excluded_drugs),
            "exclusion_active": True,
            "last_updated": "2024-09-26"
        }

# Global instance
alzheimer_filter = AlzheimerExclusionFilter()