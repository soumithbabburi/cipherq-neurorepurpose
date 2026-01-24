"""
Centralized Database Utility
All modules should use this to prevent "connection already closed" errors
"""
import psycopg2
from psycopg2 import pool
from psycopg2.extras import RealDictCursor
import os
import logging

logger = logging.getLogger(__name__)

# Global connection pool
_pool = None

def init_connection_pool():
    """Initialize connection pool (call once at startup)"""
    global _pool
    
    if _pool is not None:
        return _pool
    
    try:
        db_params = {
            "host": os.getenv("DB_HOST", "localhost"),
            "port": int(os.getenv("DB_PORT", 5432)),
            "database": os.getenv("DB_NAME", "postgres"),
            "user": os.getenv("DB_USER", "postgres"),
            "password": os.getenv("DB_PASSWORD", ""),
            "sslmode": "require" if "neon.tech" in os.getenv("DB_HOST", "") else "prefer"
        }
        
        _pool = psycopg2.pool.ThreadedConnectionPool(
            minconn=1,
            maxconn=20,
            **db_params
        )
        
        logger.info(f"✅ Database pool initialized: {db_params['database']}@{db_params['host']}")
        return _pool
        
    except Exception as e:
        logger.error(f"Failed to create connection pool: {e}")
        return None

def get_connection():
    """Get connection from pool"""
    global _pool
    
    if _pool is None:
        _pool = init_connection_pool()
    
    if _pool is None:
        raise Exception("No database connection pool available")
    
    return _pool.getconn()

def return_connection(conn):
    """Return connection to pool"""
    global _pool
    if _pool and conn:
        _pool.putconn(conn)

def execute_query(sql, params=None, fetch=True):
    """
    Execute query with automatic connection management.
    Use this in ALL modules!
    
    Example:
        results = execute_query("SELECT * FROM drugs WHERE name = %s", ('Metformin',))
    """
    conn = None
    try:
        conn = get_connection()
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        
        cursor.execute(sql, params)
        
        if fetch:
            results = cursor.fetchall()
            cursor.close()
            return [dict(row) for row in results]
        else:
            conn.commit()
            cursor.close()
            return None
            
    except Exception as e:
        logger.error(f"Query failed: {e}")
        if conn:
            conn.rollback()
        return [] if fetch else None
    finally:
        if conn:
            return_connection(conn)

def get_drug_smiles(drug_name):
    """Helper: Get SMILES for a drug"""
    results = execute_query(
        "SELECT smiles FROM drugs WHERE name = %s LIMIT 1",
        (drug_name,)
    )
    return results[0]['smiles'] if results else None

def get_protein_by_gene(gene_symbol):
    """Helper: Get protein by gene symbol"""
    results = execute_query(
        "SELECT * FROM proteins WHERE gene_symbol = %s LIMIT 1",
        (gene_symbol,)
    )
    return results[0] if results else None

# Initialize pool when module is imported
init_connection_pool()
