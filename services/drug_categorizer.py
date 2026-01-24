"""
Fixed Drug Categorizer Service
Returns ALL drugs in category (not limited to 7!)
"""
import psycopg2
import os
import logging
from typing import List, Dict

logger = logging.getLogger(__name__)

class DrugCategorizer:
    """Categorize and retrieve drugs from database - SHOWS ALL DRUGS"""
    
    def __init__(self):
        self.conn = psycopg2.connect(
            host=os.getenv("DB_HOST", "localhost"),
            database=os.getenv("DB_NAME", "cipherq_repurpose"),
            user=os.getenv("DB_USER", "babburisoumith"),
            password=os.getenv("DB_PASSWORD", "")
        )
        logger.info("Drug categorizer initialized")
    
    def get_drugs_by_category(self, category: str, limit: int = 100) -> List[Dict]:
        """
        Get ALL drugs in category (increased limit to 100!)
        
        Args:
            category: Category name (e.g., "Diabetes", "Cardiovascular")
            limit: Maximum number of drugs (default 100, was 7!)
        """
        cursor = self.conn.cursor()
        
        try:
            # Query with LARGE limit to show all drugs
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
            
            cursor.execute(sql, (f'%{category}%', limit))
            rows = cursor.fetchall()
            
            drugs = []
            for row in rows:
                drugs.append({
                    'name': row[0],
                    'category': row[1],
                    'class': row[2],
                    'mechanism': row[3],
                    'smiles': row[4],
                    'qed_score': float(row[5]) if row[5] else 0.0,
                    'fda_status': row[6]
                })
            
            cursor.close()
            logger.info(f"✅ Retrieved {len(drugs)} drugs for category {category}")
            return drugs
            
        except Exception as e:
            logger.error(f"Drug categorizer error: {e}")
            cursor.close()
            return []
    
    def get_all_categories(self) -> List[str]:
        """Get all unique categories"""
        cursor = self.conn.cursor()
        
        cursor.execute("""
            SELECT DISTINCT therapeutic_category 
            FROM drugs 
            WHERE therapeutic_category IS NOT NULL
            ORDER BY therapeutic_category
        """)
        
        categories = [row[0] for row in cursor.fetchall()]
        cursor.close()
        return categories
    
    def get_random_drugs(self, limit: int = 10) -> List[Dict]:
        """Get random drugs"""
        cursor = self.conn.cursor()
        
        cursor.execute("""
            SELECT name, therapeutic_category, drug_class, smiles
            FROM drugs
            ORDER BY RANDOM()
            LIMIT %s
        """, (limit,))
        
        rows = cursor.fetchall()
        drugs = [{'name': r[0], 'category': r[1], 'class': r[2], 'smiles': r[3]} for r in rows]
        
        cursor.close()
        return drugs


# Singleton instance
_categorizer = None

def get_drug_categorizer():
    """Get singleton DrugCategorizer instance"""
    global _categorizer
    if _categorizer is None:
        _categorizer = DrugCategorizer()
    return _categorizer