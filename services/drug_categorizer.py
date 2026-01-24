"""
Fixed Drug Categorizer Service
Uses centralized database_utils (no more connection errors!)
"""
import logging
from typing import List, Dict
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

logger = logging.getLogger(__name__)

class DrugCategorizer:
    """Categorize and retrieve drugs - uses database_utils for connections"""
    
    def __init__(self):
        # Don't create connection here - use database_utils!
        logger.info("Drug categorizer initialized (using database_utils)")
    
    def get_drugs_by_category(self, category: str, limit: int = 100) -> List[Dict]:
        """Get ALL drugs in category"""
        
        try:
            # Use centralized database connection
            from database_utils import execute_query
            
            sql = """
                SELECT 
                    name,
                    therapeutic_category,
                    drug_class,
                    mechanism_of_action,
                    smiles,
                    qed_score,
                    fda_status
                FROM drugs
                WHERE therapeutic_category ILIKE %s
                ORDER BY qed_score DESC NULLS LAST
                LIMIT %s
            """
            
            results = execute_query(sql, (f'%{category}%', limit))
            
            drugs = []
            for row in results:
                drugs.append({
                    'name': row.get('name'),
                    'category': row.get('therapeutic_category'),
                    'class': row.get('drug_class'),
                    'mechanism': row.get('mechanism_of_action'),
                    'smiles': row.get('smiles'),
                    'qed_score': float(row.get('qed_score', 0)) if row.get('qed_score') else 0.0,
                    'fda_status': row.get('fda_status')
                })
            
            logger.info(f"✅ Retrieved {len(drugs)} drugs for category {category}")
            return drugs
            
        except Exception as e:
            logger.error(f"Drug categorizer error: {e}")
            return []
    
    def get_all_categories(self) -> List[str]:
        """Get all unique categories"""
        try:
            from database_utils import execute_query
            
            results = execute_query("""
                SELECT DISTINCT therapeutic_category 
                FROM drugs 
                WHERE therapeutic_category IS NOT NULL
                ORDER BY therapeutic_category
            """)
            
            return [row['therapeutic_category'] for row in results]
        except Exception as e:
            logger.error(f"Get categories error: {e}")
            return []
    
    def get_random_drugs(self, limit: int = 10) -> List[Dict]:
        """Get random drugs"""
        try:
            from database_utils import execute_query
            
            results = execute_query("""
                SELECT name, therapeutic_category, drug_class, smiles
                FROM drugs
                ORDER BY RANDOM()
                LIMIT %s
            """, (limit,))
            
            return [{'name': r['name'], 'category': r['therapeutic_category'], 
                    'class': r['drug_class'], 'smiles': r['smiles']} for r in results]
        except Exception as e:
            logger.error(f"Random drugs error: {e}")
            return []


# Singleton instance
_categorizer = None

def get_drug_categorizer():
    """Get singleton DrugCategorizer instance"""
    global _categorizer
    if _categorizer is None:
        _categorizer = DrugCategorizer()
    return _categorizer
