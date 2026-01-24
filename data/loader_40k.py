"""
Load 100k+ datasets with REAL drug data
Shows: Real names, Valid SMILES, Indications, Real targets
"""

import json
import logging
import os
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# Import centralized database queries
try:
    from database_queries import (
        get_db_connection, get_all_drugs, get_drug_by_name,
        get_drugs_by_category, search_drugs, get_drug_targets
    )
    DATABASE_QUERIES_AVAILABLE = True
except ImportError:
    DATABASE_QUERIES_AVAILABLE = False
    logger.warning("database_queries not available, using fallback")


class Data40kLoader:
    """Load and cache 100k+ biomedical datasets - uses PostgreSQL"""
    
    def __init__(self):
        self.data_dir = Path("data")
        self.assets_dir = Path("attached_assets")
        self._drugs = None
        self._genes = None
        self._proteins = None
        self._pathways = None
        self._interactions = None
        self._use_database = True
    
    def _load_drugs_from_database(self) -> List[Dict]:
        """Load REAL drugs from PostgreSQL with indications"""
        if not DATABASE_QUERIES_AVAILABLE:
            return None
        
        try:
            # Load ALL drugs (100k with REAL data)
            drugs_data = get_all_drugs(limit=100000)
            
            if not drugs_data:
                return None
            
            drugs = []
            for row in drugs_data:
                drug_name = row['name']
                
                # Get real protein targets
                targets = get_drug_targets(drug_name, limit=5)
                
                # Format target string with REAL proteins
                if targets:
                    target_genes = [t['gene_symbol'] for t in targets[:3]]
                    target_str = ', '.join(target_genes)
                    if len(targets) > 3:
                        target_str += f" (+{len(targets)-3})"
                else:
                    target_str = 'Unknown'
                
                # Calculate stats
                avg_confidence = 0
                avg_binding = None
                if targets:
                    avg_confidence = sum(t.get('confidence_score', 0) for t in targets) / len(targets)
                    bindings = [t.get('binding_affinity') for t in targets if t.get('binding_affinity')]
                    if bindings:
                        avg_binding = sum(bindings) / len(bindings)
                
                # Format drug with ALL data
                drugs.append({
                    'id': row['id'],
                    'name': drug_name,
                    'class': row.get('drug_class', 'Unknown'),
                    'therapeutic_category': row.get('therapeutic_category', 'General'),
                    'indication': row.get('original_indication', 'Multiple indications'),  # NEW!
                    'target': target_str,  # REAL targets
                    'mechanism': row.get('mechanism_of_action') or 'Under investigation',
                    'smiles': row.get('smiles', ''),
                    'source': 'CipherQ Database',
                    'status': row.get('fda_status', 'Unknown'),
                    'qed_score': row.get('qed_score'),
                    'target_count': len(targets),
                    'confidence': round(avg_confidence, 3) if avg_confidence else 0,
                    'binding_affinity': round(avg_binding, 2) if avg_binding else None,
                    'top_targets': [t['gene_symbol'] for t in targets[:5]],
                    # Additional useful fields
                    'molecular_weight': row.get('molecular_weight'),
                    'log_p': row.get('log_p'),
                    'oral_bioavailability': row.get('oral_bioavailability')
                })
            
            logger.info(f"Loaded {len(drugs)} REAL drugs from PostgreSQL database")
            logger.info(f"  With: Valid SMILES, Indications, Real targets")
            return drugs
        except Exception as e:
            logger.warning(f"Error loading drugs from database: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None
    
    def _load_drugs_from_json(self) -> List[Dict]:
        """Fallback: Load drugs from JSON file"""
        filepath = self.data_dir / "drugs_40k.json"
        try:
            with open(filepath, 'r') as f:
                drugs = json.load(f)
            logger.info(f"Loaded {len(drugs)} drugs from drugs_40k.json (fallback)")
            return drugs
        except Exception as e:
            logger.error(f"Error loading drugs_40k.json: {e}")
            return []
    
    @property
    def drugs(self) -> List[Dict]:
        """Load drugs - database first, fallback to JSON"""
        if self._drugs is None:
            if self._use_database:
                self._drugs = self._load_drugs_from_database()
            
            if self._drugs is None:
                self._drugs = self._load_drugs_from_json()
        
        return self._drugs
    
    @property
    def genes(self) -> List[Dict]:
        """Load genes from assets"""
        if self._genes is None:
            filepath = self.assets_dir / "genes_40k_1763140250778.json"
            try:
                with open(filepath, 'r') as f:
                    self._genes = json.load(f)
                logger.info(f"Loaded {len(self._genes)} genes from genes dataset")
            except Exception as e:
                logger.error(f"Error loading genes: {e}")
                self._genes = []
        return self._genes
    
    @property
    def proteins(self) -> List[Dict]:
        """Load proteins"""
        if self._proteins is None:
            filepath_40k = self.assets_dir / "proteins_40k_1763140250779.json"
            filepath_500 = self.data_dir / "proteins_500.json"
            
            try:
                if filepath_40k.exists():
                    with open(filepath_40k, 'r') as f:
                        self._proteins = json.load(f)
                    logger.info(f"Loaded {len(self._proteins)} proteins from 40k dataset")
                elif filepath_500.exists():
                    with open(filepath_500, 'r') as f:
                        self._proteins = json.load(f)
                    logger.info(f"Loaded {len(self._proteins)} proteins from 500 dataset")
                else:
                    logger.error("No protein data files found")
                    self._proteins = []
            except Exception as e:
                logger.error(f"Error loading proteins: {e}")
                self._proteins = []
        return self._proteins
    
    @property
    def pathways(self) -> List[Dict]:
        """Load pathways"""
        if self._pathways is None:
            filepath_40k = self.assets_dir / "pathways_40k_1763140250779.json"
            filepath_500 = self.data_dir / "pathways_500.json"
            
            try:
                if filepath_40k.exists():
                    with open(filepath_40k, 'r') as f:
                        self._pathways = json.load(f)
                    logger.info(f"Loaded {len(self._pathways)} pathways from 40k dataset")
                elif filepath_500.exists():
                    with open(filepath_500, 'r') as f:
                        self._pathways = json.load(f)
                    logger.info(f"Loaded {len(self._pathways)} pathways from 500 dataset")
                else:
                    logger.error("No pathway data files found")
                    self._pathways = []
            except Exception as e:
                logger.error(f"Error loading pathways: {e}")
                self._pathways = []
        return self._pathways
    
    @property
    def interactions(self) -> List[Dict]:
        """Load interactions"""
        if self._interactions is None:
            filepath = self.assets_dir / "interactions_40k_1763140250779.json"
            try:
                with open(filepath, 'r') as f:
                    self._interactions = json.load(f)
                logger.info(f"Loaded {len(self._interactions)} interactions")
            except Exception as e:
                logger.error(f"Error loading interactions: {e}")
                self._interactions = []
        return self._interactions
    
    def get_drug_by_name(self, drug_name: str) -> Optional[Dict]:
        """Get drug by name with all details"""
        if self._use_database and DATABASE_QUERIES_AVAILABLE:
            try:
                return get_drug_by_name(drug_name)
            except Exception as e:
                logger.warning(f"Database query failed: {e}")
        
        # Fallback
        for drug in self.drugs:
            if drug.get('name', '').lower() == drug_name.lower():
                return drug
        return None
    
    def get_all_drug_names(self) -> List[str]:
        """Get all drug names"""
        return [drug.get('name') for drug in self.drugs if drug.get('name')]
    
    def get_drugs_by_category(self, category: str) -> List[Dict]:
        """Get drugs by category"""
        if self._use_database and DATABASE_QUERIES_AVAILABLE:
            try:
                return get_drugs_by_category(category, limit=1000)
            except:
                pass
        
        return [d for d in self.drugs if d.get('therapeutic_category') == category]
    
    def search_drugs(self, query: str, limit: int = 50) -> List[Dict]:
        """Search drugs"""
        if self._use_database and DATABASE_QUERIES_AVAILABLE:
            try:
                return search_drugs(query, limit=limit)
            except:
                pass
        
        query_lower = query.lower()
        matches = [d for d in self.drugs if query_lower in d.get('name', '').lower()]
        return matches[:limit]
    
    def get_interactions_for_drug(self, drug_name: str) -> List[Dict]:
        """Get interactions for drug - uses database"""
        if DATABASE_QUERIES_AVAILABLE:
            try:
                return get_drug_targets(drug_name, limit=20)
            except:
                pass
        
        # Fallback to JSON
        drug_upper = drug_name.upper()
        return [i for i in self.interactions if i.get('drug_name', '').upper() == drug_upper]


data_40k = Data40kLoader()