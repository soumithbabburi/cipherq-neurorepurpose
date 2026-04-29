"""
NeuroRepurpose Database Schema
Manages the neurorepurpose PostgreSQL database backed by ChEMBL 33 + HetioNet data.
"""

import psycopg2
import psycopg2.pool
import logging
import os
from contextlib import contextmanager

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

DB_CONFIG = {
    "host": os.getenv("DB_HOST", "localhost"),
    "port": int(os.getenv("DB_PORT", "5434")),
    "user": os.getenv("DB_USER", "babburisoumith"),
    "password": os.getenv("DB_PASSWORD", ""),
    "dbname": "neurorepurpose",
}

_pool = None


def get_pool():
    global _pool
    if _pool is None:
        try:
            _pool = psycopg2.pool.ThreadedConnectionPool(
                minconn=1,
                maxconn=10,
                **DB_CONFIG,
            )
        except Exception as e:
            logger.error(f"Failed to create connection pool: {e}")
            raise
    return _pool


@contextmanager
def get_connection():
    pool = get_pool()
    conn = pool.getconn()
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        pool.putconn(conn)


def _create_database():
    """Create the neurorepurpose database if it doesn't exist."""
    admin_cfg = {**DB_CONFIG, "dbname": "postgres"}
    conn = psycopg2.connect(**admin_cfg)
    conn.autocommit = True
    cur = conn.cursor()
    cur.execute("SELECT 1 FROM pg_database WHERE datname = 'neurorepurpose'")
    if not cur.fetchone():
        cur.execute("CREATE DATABASE neurorepurpose")
        logger.info("Created database: neurorepurpose")
    else:
        logger.info("Database neurorepurpose already exists")
    conn.close()


DDL = """
-- ─── Compounds ────────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS compounds (
    id              SERIAL PRIMARY KEY,
    chembl_id       VARCHAR(20) UNIQUE,
    name            TEXT NOT NULL,
    smiles          TEXT,
    max_phase       REAL DEFAULT 0,
    therapeutic_flag BOOLEAN DEFAULT FALSE,
    dosed_ingredient BOOLEAN DEFAULT FALSE,
    source          VARCHAR(50) DEFAULT 'chembl',
    created_at      TIMESTAMP DEFAULT NOW()
);
CREATE INDEX IF NOT EXISTS idx_compounds_name ON compounds (LOWER(name));
CREATE INDEX IF NOT EXISTS idx_compounds_chembl ON compounds (chembl_id);

-- ─── Compound Properties ──────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS compound_properties (
    compound_id     INTEGER REFERENCES compounds(id) ON DELETE CASCADE,
    mw              REAL,
    alogp           REAL,
    psa             REAL,
    hba             INTEGER,
    hbd             INTEGER,
    rtb             INTEGER,
    ro5_violations  INTEGER DEFAULT 0,
    qed             REAL,
    PRIMARY KEY (compound_id)
);

-- ─── Targets ─────────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS targets (
    id              SERIAL PRIMARY KEY,
    chembl_tid      VARCHAR(20) UNIQUE,
    name            TEXT NOT NULL,
    target_type     VARCHAR(100),
    gene_symbol     VARCHAR(50),
    organism        VARCHAR(100) DEFAULT 'Homo sapiens'
);
CREATE INDEX IF NOT EXISTS idx_targets_name ON targets (LOWER(name));
CREATE INDEX IF NOT EXISTS idx_targets_gene ON targets (gene_symbol);

-- ─── Mechanisms ───────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS mechanisms (
    id              SERIAL PRIMARY KEY,
    compound_id     INTEGER REFERENCES compounds(id) ON DELETE CASCADE,
    target_id       INTEGER REFERENCES targets(id) ON DELETE CASCADE,
    mechanism       TEXT,
    action_type     VARCHAR(100),
    confidence      REAL DEFAULT 0.5
);
CREATE INDEX IF NOT EXISTS idx_mech_compound ON mechanisms (compound_id);
CREATE INDEX IF NOT EXISTS idx_mech_target   ON mechanisms (target_id);

-- ─── Disease Indications ──────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS indications (
    id              SERIAL PRIMARY KEY,
    compound_id     INTEGER REFERENCES compounds(id) ON DELETE CASCADE,
    disease         TEXT NOT NULL,
    efo_id          VARCHAR(50),
    max_phase       REAL,
    source          VARCHAR(50) DEFAULT 'chembl'
);
CREATE INDEX IF NOT EXISTS idx_ind_compound ON indications (compound_id);
CREATE INDEX IF NOT EXISTS idx_ind_disease  ON indications (LOWER(disease));

-- ─── HetioNet Nodes ───────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS hetionet_nodes (
    id      TEXT PRIMARY KEY,
    name    TEXT NOT NULL,
    kind    VARCHAR(100) NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_hn_kind ON hetionet_nodes (kind);
CREATE INDEX IF NOT EXISTS idx_hn_name ON hetionet_nodes (LOWER(name));

-- ─── HetioNet Edges ───────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS hetionet_edges (
    id          SERIAL PRIMARY KEY,
    source_id   TEXT REFERENCES hetionet_nodes(id) ON DELETE CASCADE,
    target_id   TEXT REFERENCES hetionet_nodes(id) ON DELETE CASCADE,
    metaedge    VARCHAR(20) NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_he_source    ON hetionet_edges (source_id);
CREATE INDEX IF NOT EXISTS idx_he_target    ON hetionet_edges (target_id);
CREATE INDEX IF NOT EXISTS idx_he_metaedge  ON hetionet_edges (metaedge);
"""


def initialize_schema():
    """Create database and all tables. Safe to call multiple times."""
    _create_database()
    with get_connection() as conn:
        cur = conn.cursor()
        cur.execute(DDL)
    logger.info("Schema initialised successfully")


if __name__ == "__main__":
    initialize_schema()
    print("neurorepurpose schema ready.")
