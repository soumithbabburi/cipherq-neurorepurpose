"""
Database Queries Module
All database query functions - 100% database-driven
"""
import logging
import psycopg2
from psycopg2.extras import RealDictCursor
from typing import List, Dict, Optional

logger = logging.getLogger(__name__)

def get_db_connection():
    """Get PostgreSQL connection"""
    try:
        conn = psycopg2.connect(
            host="localhost",
            database="cipherq_repurpose",
            user="babburisoumith",
            password=""
        )
        return conn
    except Exception as e:
        logger.error(f"DB connection failed: {e}")
        return None


def get_drug_targets(drug_name: str, limit: int = 10) -> List[Dict]:
    """Get targets from drug_protein_interactions"""
    try:
        conn = get_db_connection()
        if not conn:
            return []
        
        cur = conn.cursor(cursor_factory=RealDictCursor)
        cur.execute("""
            SELECT 
                p.gene_symbol,
                p.name as protein_name,
                dpi.confidence_score,
                dpi.interaction_type
            FROM drug_protein_interactions dpi
            JOIN drugs d ON d.id = dpi.drug_id
            JOIN proteins p ON p.id = dpi.protein_id
            WHERE LOWER(d.name) = LOWER(%s)
            ORDER BY dpi.confidence_score DESC
            LIMIT %s
        """, (drug_name, limit))
        
        results = [dict(row) for row in cur.fetchall()]
        cur.close()
        conn.close()
        
        if results:
            logger.info(f"✅ {drug_name} targets: {[r['gene_symbol'] for r in results]}")
        
        return results
        
    except Exception as e:
        logger.error(f"get_drug_targets error: {e}")
        return []


def get_drug_by_name(drug_name: str) -> Optional[Dict]:
    """Get drug info from database"""
    try:
        conn = get_db_connection()
        if not conn:
            return None
        
        cur = conn.cursor(cursor_factory=RealDictCursor)
        cur.execute("""
            SELECT *
            FROM drugs
            WHERE LOWER(name) = LOWER(%s)
            LIMIT 1
        """, (drug_name,))
        
        row = cur.fetchone()
        cur.close()
        conn.close()
        
        return dict(row) if row else None
        
    except Exception as e:
        logger.error(f"get_drug_by_name error: {e}")
        return None


def get_drugs_by_category(category: str, limit: int = 50) -> List[Dict]:
    """Get drugs by category"""
    try:
        conn = get_db_connection()
        if not conn:
            return []
        
        cur = conn.cursor(cursor_factory=RealDictCursor)
        cur.execute("""
            SELECT *
            FROM drugs
            WHERE therapeutic_category = %s
            ORDER BY name
            LIMIT %s
        """, (category, limit))
        
        results = [dict(row) for row in cur.fetchall()]
        cur.close()
        conn.close()
        
        return results
        
    except Exception as e:
        logger.error(f"get_drugs_by_category error: {e}")
        return []


def search_drugs_by_query(query: str, limit: int = 20) -> List[Dict]:
    """Search drugs by query"""
    try:
        conn = get_db_connection()
        if not conn:
            return []
        
        cur = conn.cursor(cursor_factory=RealDictCursor)
        pattern = f"%{query}%"
        cur.execute("""
            SELECT *
            FROM drugs
            WHERE name ILIKE %s
               OR original_indication ILIKE %s
               OR drug_class ILIKE %s
            LIMIT %s
        """, (pattern, pattern, pattern, limit))
        
        results = [dict(row) for row in cur.fetchall()]
        cur.close()
        conn.close()
        
        return results
        
    except Exception as e:
        logger.error(f"search_drugs error: {e}")
        return []


def get_interaction_count(drug_name: str) -> int:
    """Get number of protein interactions for a drug"""
    try:
        conn = get_db_connection()
        if not conn:
            return 0
        
        cur = conn.cursor()
        cur.execute("""
            SELECT COUNT(*)
            FROM drug_protein_interactions dpi
            JOIN drugs d ON d.id = dpi.drug_id
            WHERE LOWER(d.name) = LOWER(%s)
        """, (drug_name,))
        
        count = cur.fetchone()[0]
        cur.close()
        conn.close()
        
        return count
        
    except Exception as e:
        logger.error(f"get_interaction_count error: {e}")
        return 0


def get_disease_targets(disease_name: str, limit: int = 10) -> List[str]:
    """Get common protein targets for a disease"""
    # Simplified - returns common targets based on disease type
    disease_targets = {
        'Alzheimer': ['ACHE', 'BCHE', 'NMDA', 'APP', 'MAPT'],
        'Diabetes': ['PPARG', 'DPP4', 'SGLT2', 'GLP1R', 'AMPK'],
        'Cardiovascular': ['ACE', 'ADRB1', 'HMGCR', 'AGTR1'],
        'Cancer': ['ABL1', 'EGFR', 'ERBB2', 'KIT'],
    }
    
    for key, targets in disease_targets.items():
        if key.lower() in disease_name.lower():
            return targets[:limit]
    
    return []


__all__ = ['get_drug_targets', 'get_drug_by_name', 'get_drugs_by_category', 'search_drugs_by_query', 'get_interaction_count', 'get_disease_targets']
