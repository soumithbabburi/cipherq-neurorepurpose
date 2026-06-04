"""
NeuroRepurpose Data Importer
Imports from: ChEMBL 33 (local PostgreSQL) + HetioNet TSV + existing JSON files.
Run standalone: python3 database/importer.py
"""

import json
import hashlib
import logging
import os
import re
import sys
from pathlib import Path

import psycopg2
import psycopg2.extras

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parent.parent

def _conn_config(prefix: str, default_name: str) -> dict:
    return {
        "host": os.getenv(f"{prefix}_HOST", "localhost"),
        "port": int(os.getenv(f"{prefix}_PORT", "5433")),
        "user": os.getenv(f"{prefix}_USER", "babburisoumith"),
        "password": os.getenv(f"{prefix}_PASSWORD", ""),
        "dbname": os.getenv(f"{prefix}_NAME", default_name),
    }


# ChEMBL connection (read-only source)
CHEMBL_CONN = _conn_config("CHEMBL_DB", "chembl_33")

# Destination connection
NEURO_CONN = {
    "host": os.getenv("DB_HOST", CHEMBL_CONN["host"]),
    "port": int(os.getenv("DB_PORT", str(CHEMBL_CONN["port"]))),
    "user": os.getenv("DB_USER", CHEMBL_CONN["user"]),
    "password": os.getenv("DB_PASSWORD", CHEMBL_CONN["password"]),
    "dbname": os.getenv("DB_NAME", "neurorepurpose"),
}

NEURO_DISEASES = [
    "%alzheimer%", "%parkinson%", "%multiple sclerosis%", "%epilep%",
    "%schizophrenia%", "%depression%", "%anxiety%", "%autism%",
    "%huntington%", "%dementia%", "%amyotrophic%", "%stroke%",
    "%bipolar%", "%psychosis%", "%neurolog%", "%neuropath%",
    "%migraine%", "%insomnia%", "%mania%", "%obsessive%",
]


def import_chembl_neuro_compounds():
    """Pull neuro-relevant compounds + targets + mechanisms + indications from ChEMBL 33."""
    logger.info("Importing ChEMBL neuro compounds…")

    like_clauses = " OR ".join(f"LOWER(di.mesh_heading) LIKE '{d}'" for d in NEURO_DISEASES)

    sql = f"""
        SELECT DISTINCT
            md.molregno,
            md.chembl_id,
            COALESCE(md.pref_name, md.chembl_id) AS name,
            COALESCE(md.max_phase, 0) AS max_phase,
            COALESCE(md.therapeutic_flag, 0)::smallint AS therapeutic_flag,
            COALESCE(md.dosed_ingredient, 0)::smallint AS dosed_ingredient,
            cs.canonical_smiles,
            cp.mw_freebase,
            cp.alogp,
            cp.psa,
            cp.hba,
            cp.hbd,
            cp.rtb,
            COALESCE(cp.num_ro5_violations, 0) AS ro5_violations
        FROM molecule_dictionary md
        JOIN drug_indication di ON di.molregno = md.molregno
        LEFT JOIN compound_structures cs ON cs.molregno = md.molregno
        LEFT JOIN compound_properties cp ON cp.molregno = md.molregno
        WHERE {like_clauses}
    """

    src = psycopg2.connect(**CHEMBL_CONN)
    dst = psycopg2.connect(**NEURO_CONN)
    src_cur = src.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    dst_cur = dst.cursor()

    src_cur.execute(sql)
    rows = src_cur.fetchall()
    logger.info(f"  ChEMBL returned {len(rows):,} neuro compounds")

    compound_map = {}  # molregno → internal id

    for r in rows:
        dst_cur.execute(
            """
            INSERT INTO compounds (chembl_id, name, smiles, max_phase, therapeutic_flag, dosed_ingredient, source)
            VALUES (%s, %s, %s, %s, %s, %s, 'chembl')
            ON CONFLICT (chembl_id) DO UPDATE
                SET name = EXCLUDED.name,
                    smiles = COALESCE(EXCLUDED.smiles, compounds.smiles),
                    max_phase = GREATEST(EXCLUDED.max_phase, compounds.max_phase)
            RETURNING id
            """,
            (r["chembl_id"], r["name"], r["canonical_smiles"], float(r["max_phase"] or 0),
             bool(r["therapeutic_flag"]), bool(r["dosed_ingredient"])),
        )
        cid = dst_cur.fetchone()[0]
        compound_map[r["molregno"]] = cid

        if any(r[k] is not None for k in ["mw_freebase", "alogp", "psa"]):
            dst_cur.execute(
                """
                INSERT INTO compound_properties (compound_id, mw, alogp, psa, hba, hbd, rtb, ro5_violations)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (compound_id) DO UPDATE
                    SET mw=EXCLUDED.mw, alogp=EXCLUDED.alogp, psa=EXCLUDED.psa,
                        hba=EXCLUDED.hba, hbd=EXCLUDED.hbd, rtb=EXCLUDED.rtb,
                        ro5_violations=EXCLUDED.ro5_violations
                """,
                (cid, r["mw_freebase"], r["alogp"], r["psa"],
                 r["hba"], r["hbd"], r["rtb"], r["ro5_violations"]),
            )

    dst.commit()
    logger.info(f"  Upserted {len(compound_map):,} compounds")

    # ── Mechanisms + targets ─────────────────────────────────────────────────
    if compound_map:
        molregnos = list(compound_map.keys())
        placeholders = ",".join(["%s"] * len(molregnos))
        src_cur.execute(
            f"""
            SELECT dm.molregno, dm.mechanism_of_action, dm.action_type,
                   td.chembl_id AS target_chembl, td.pref_name AS target_name,
                   td.target_type, td.organism
            FROM drug_mechanism dm
            JOIN target_dictionary td ON td.tid = dm.tid
            WHERE dm.molregno IN ({placeholders})
            """,
            molregnos,
        )
        mechs = src_cur.fetchall()
        logger.info(f"  Got {len(mechs):,} mechanisms")
        target_cache = {}

        for m in mechs:
            tc = m["target_chembl"]
            if tc not in target_cache:
                dst_cur.execute(
                    """
                    INSERT INTO targets (chembl_tid, name, target_type, organism)
                    VALUES (%s, %s, %s, %s)
                    ON CONFLICT (chembl_tid) DO UPDATE SET name=EXCLUDED.name
                    RETURNING id
                    """,
                    (tc, m["target_name"], m["target_type"], m["organism"]),
                )
                target_cache[tc] = dst_cur.fetchone()[0]

            cid = compound_map.get(m["molregno"])
            if cid:
                dst_cur.execute(
                    """
                    INSERT INTO mechanisms (compound_id, target_id, mechanism, action_type, confidence)
                    VALUES (%s, %s, %s, %s, 0.9)
                    ON CONFLICT DO NOTHING
                    """,
                    (cid, target_cache[tc], m["mechanism_of_action"], m["action_type"]),
                )

        dst.commit()

    # ── Indications ──────────────────────────────────────────────────────────
    if compound_map:
        src_cur.execute(
            f"""
            SELECT di.molregno, di.mesh_heading, di.efo_id, di.max_phase_for_ind
            FROM drug_indication di
            WHERE di.molregno IN ({placeholders})
            """,
            molregnos,
        )
        inds = src_cur.fetchall()
        logger.info(f"  Got {len(inds):,} indications")
        for ind in inds:
            cid = compound_map.get(ind["molregno"])
            if cid:
                dst_cur.execute(
                    """
                    INSERT INTO indications (compound_id, disease, efo_id, max_phase, source)
                    VALUES (%s, %s, %s, %s, 'chembl')
                    ON CONFLICT DO NOTHING
                    """,
                    (cid, ind["mesh_heading"], ind["efo_id"], ind["max_phase_for_ind"]),
                )
        dst.commit()

    src.close()
    dst.close()
    logger.info("ChEMBL import complete.")


def import_hetionet_data():
    """Import HetioNet nodes and edges TSV files."""
    nodes_path = REPO_ROOT / "data" / "hetionet" / "nodes.tsv"
    edges_path = REPO_ROOT / "data" / "hetionet" / "edges.tsv"

    if not nodes_path.exists():
        logger.warning("HetioNet nodes.tsv not found — skipping")
        return

    conn = psycopg2.connect(**NEURO_CONN)
    cur = conn.cursor()

    # Nodes
    logger.info("Importing HetioNet nodes…")
    nodes = []
    with open(nodes_path, encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                nodes.append((parts[0], parts[1], parts[2]))

    psycopg2.extras.execute_values(
        cur,
        "INSERT INTO hetionet_nodes (id, name, kind) VALUES %s ON CONFLICT (id) DO NOTHING",
        nodes,
        page_size=1000,
    )
    conn.commit()
    logger.info(f"  Imported {len(nodes):,} HetioNet nodes")

    # Edges
    if not edges_path.exists():
        logger.warning("HetioNet edges.tsv not found — skipping edges")
        conn.close()
        return

    # Verify it's real TSV not HTML
    with open(edges_path, "rb") as f:
        header_bytes = f.read(100)
    if b"<html" in header_bytes.lower() or b"404" in header_bytes:
        logger.warning("HetioNet edges.tsv appears to be HTML/error — skipping")
        conn.close()
        return

    logger.info("Importing HetioNet edges…")
    edges = []
    with open(edges_path, encoding="utf-8") as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                edges.append((parts[0], parts[2], parts[1]))  # source, target, metaedge

    batch_size = 5000
    total = 0
    for i in range(0, len(edges), batch_size):
        batch = edges[i : i + batch_size]
        try:
            psycopg2.extras.execute_values(
                cur,
                "INSERT INTO hetionet_edges (source_id, target_id, metaedge) VALUES %s ON CONFLICT DO NOTHING",
                batch,
            )
            conn.commit()
            total += len(batch)
        except Exception as e:
            conn.rollback()
            logger.warning(f"Edge batch {i}-{i+batch_size} failed: {e}")

    logger.info(f"  Imported {total:,} HetioNet edges")
    conn.close()


def _phase_to_num(phase) -> float:
    text = str(phase or "").lower()
    if "launched" in text or "approved" in text or "phase 4" in text:
        return 4.0
    if "phase 3" in text:
        return 3.0
    if "phase 2" in text:
        return 2.0
    if "phase 1" in text:
        return 1.0
    return 0.0


def _norm_name(name: str) -> str:
    return re.sub(r"\s+", " ", str(name or "").strip()).lower()


def import_chembl_snapshot_files():
    """Import ChEMBL-derived snapshot JSON files bundled with the repo."""
    logger.info("Importing bundled ChEMBL snapshot files…")
    snapshot_dir = REPO_ROOT / "data_small_molecules"
    files = [
        snapshot_dir / "chembl_launched.json",
        snapshot_dir / "chembl_phase3.json",
        snapshot_dir / "chembl_phase2.json",
        snapshot_dir / "chembl_phase1.json",
    ]

    conn = psycopg2.connect(**NEURO_CONN)
    cur = conn.cursor()
    imported = 0
    for path in files:
        if not path.exists():
            logger.warning(f"  Missing {path.name}; skipping")
            continue
        with open(path, encoding="utf-8") as f:
            rows = json.load(f)

        for row in rows:
            name = row.get("name")
            if not name:
                continue
            phase = _phase_to_num(row.get("clinical_phase"))
            chembl_id = row.get("chembl_id") or f"C15_{hashlib.sha1(_norm_name(name).encode('utf-8')).hexdigest()[:16]}"
            cur.execute(
                """
                INSERT INTO compounds (chembl_id, name, smiles, max_phase, therapeutic_flag, dosed_ingredient, source)
                VALUES (%s, %s, %s, %s, %s, %s, 'chembl_snapshot')
                ON CONFLICT (chembl_id) DO UPDATE
                    SET name=EXCLUDED.name,
                        smiles=COALESCE(EXCLUDED.smiles, compounds.smiles),
                        max_phase=GREATEST(EXCLUDED.max_phase, compounds.max_phase),
                        source='chembl_snapshot'
                RETURNING id
                """,
                (chembl_id, name, row.get("smiles"), phase, phase >= 3, phase >= 4),
            )
            cid = cur.fetchone()[0]
            cur.execute(
                """
                INSERT INTO compound_properties (compound_id, mw, alogp, psa, hba, hbd, rtb, ro5_violations)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (compound_id) DO UPDATE
                    SET mw=EXCLUDED.mw, alogp=EXCLUDED.alogp, psa=EXCLUDED.psa,
                        hba=EXCLUDED.hba, hbd=EXCLUDED.hbd, rtb=EXCLUDED.rtb,
                        ro5_violations=EXCLUDED.ro5_violations
                """,
                (
                    cid,
                    row.get("mw"),
                    row.get("logp"),
                    row.get("psa"),
                    row.get("hba"),
                    row.get("hbd"),
                    row.get("rotatable_bonds"),
                    row.get("lipinski_violations") or 0,
                ),
            )
            if row.get("pchembl_value") is not None:
                cur.execute(
                    """
                    INSERT INTO compound_activities
                        (compound_id, target_id, activity_type, pchembl_value, standard_value, standard_units)
                    VALUES (%s, NULL, %s, %s, NULL, NULL)
                    """,
                    (cid, row.get("assay_type") or "activity", row.get("pchembl_value")),
                )
            imported += 1

    conn.commit()
    conn.close()
    logger.info(f"  Imported/updated {imported:,} ChEMBL snapshot compounds")


def import_json_associations_and_interactions():
    """Backfill indications and target mechanisms from repo data files."""
    logger.info("Importing data-driven indications and target interactions…")
    assoc_path = REPO_ROOT / "disease_associations.json"
    interactions_path = REPO_ROOT / "drug_interactions.json"
    associations = json.load(open(assoc_path, encoding="utf-8")) if assoc_path.exists() else {}
    interactions = json.load(open(interactions_path, encoding="utf-8")) if interactions_path.exists() else {}

    conn = psycopg2.connect(**NEURO_CONN)
    cur = conn.cursor()
    cur.execute("SELECT id, LOWER(name), max_phase FROM compounds")
    compound_map = {row[1]: (row[0], float(row[2] or 0)) for row in cur.fetchall()}
    target_cache = {}

    ind_count = 0
    for drug_name, diseases in associations.items():
        match = compound_map.get(_norm_name(drug_name))
        if not match:
            continue
        cid, phase = match
        for disease in diseases:
            cur.execute(
                """
                INSERT INTO indications (compound_id, disease, mesh_id, efo_id, max_phase, source)
                VALUES (%s, %s, NULL, NULL, %s, 'repo_associations')
                ON CONFLICT DO NOTHING
                """,
                (cid, disease, phase),
            )
            ind_count += cur.rowcount

    mech_count = 0
    for drug_name, rows in interactions.items():
        match = compound_map.get(_norm_name(drug_name))
        if not match:
            continue
        cid, _phase = match
        for item in rows:
            gene = item.get("gene_symbol")
            target_name = item.get("protein_name") or gene
            if not target_name:
                continue
            key = gene or target_name
            if key not in target_cache:
                target_id = f"RT_{hashlib.sha1(str(key).encode('utf-8')).hexdigest()[:17]}"
                cur.execute(
                    """
                    INSERT INTO targets (chembl_tid, name, target_type, gene_symbol, organism)
                    VALUES (%s, %s, 'Protein', %s, 'Homo sapiens')
                    ON CONFLICT (chembl_tid) DO UPDATE SET name=EXCLUDED.name, gene_symbol=EXCLUDED.gene_symbol
                    RETURNING id
                    """,
                    (target_id, target_name, gene),
                )
                target_cache[key] = cur.fetchone()[0]
            cur.execute(
                """
                INSERT INTO mechanisms (compound_id, target_id, mechanism, action_type, confidence)
                VALUES (%s, %s, %s, %s, %s)
                ON CONFLICT DO NOTHING
                """,
                (
                    cid,
                    target_cache[key],
                    item.get("interaction_type") or "target interaction",
                    item.get("interaction_type") or "target",
                    item.get("confidence_score") or 0.5,
                ),
            )
            mech_count += cur.rowcount

    conn.commit()
    conn.close()
    logger.info(f"  Imported {ind_count:,} indications and {mech_count:,} mechanisms")


def import_json_drugs():
    """Import drugs from existing JSON files that aren't already in the database."""
    drugs_path = REPO_ROOT / "drugs.json"
    if not drugs_path.exists():
        logger.warning("drugs.json not found — skipping JSON import")
        return

    logger.info("Importing drugs from drugs.json…")
    with open(drugs_path, encoding="utf-8") as f:
        drugs = json.load(f)

    conn = psycopg2.connect(**NEURO_CONN)
    cur = conn.cursor()

    inserted = 0
    for name, data in drugs.items():
        smiles = data.get("smiles", "")
        mw = data.get("molecular_weight")
        logp = data.get("logp", data.get("log_p"))
        max_phase = 4.0 if data.get("fda_approved") or data.get("approved") else 0.0

        cur.execute(
            """
            INSERT INTO compounds (name, smiles, max_phase, source)
            VALUES (%s, %s, %s, 'json')
            ON CONFLICT DO NOTHING
            RETURNING id
            """,
            (name.title(), smiles, max_phase),
        )
        row = cur.fetchone()
        if row:
            cid = row[0]
            if mw is not None or logp is not None:
                cur.execute(
                    """
                    INSERT INTO compound_properties (compound_id, mw, alogp)
                    VALUES (%s, %s, %s)
                    ON CONFLICT DO NOTHING
                    """,
                    (cid, mw, logp),
                )
            inserted += 1

    conn.commit()
    conn.close()
    logger.info(f"  Imported {inserted:,} new drugs from JSON")


def generate_chembl_edges():
    """
    Generate HetioNet-style edges from ChEMBL compound-target-disease data.
    Uses the nodes already imported from hetionet-v1.0-nodes.tsv as anchors
    where possible, otherwise creates new nodes.
    Metaedges used:
      CtD  Compound treats Disease
      CbG  Compound binds Gene  (from drug_mechanism targets)
      GaD  Gene associates Disease (inferred via shared compound)
    """
    logger.info("Generating knowledge graph edges from ChEMBL data…")

    conn = psycopg2.connect(**NEURO_CONN)
    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    w_cur = conn.cursor()

    # Fetch all compounds with their targets and indications
    cur.execute("""
        SELECT c.id, c.name, c.chembl_id,
               STRING_AGG(DISTINCT i.disease, '||') AS diseases,
               STRING_AGG(DISTINCT t.name, '||') AS targets,
               STRING_AGG(DISTINCT t.gene_symbol, '||') AS genes
        FROM compounds c
        LEFT JOIN indications i ON i.compound_id = c.id
        LEFT JOIN mechanisms m ON m.compound_id = c.id
        LEFT JOIN targets t ON t.id = m.target_id
        WHERE c.source IN ('chembl', 'chembl_snapshot', 'json')
        GROUP BY c.id, c.name, c.chembl_id
    """)
    rows = cur.fetchall()
    logger.info(f"  Processing {len(rows):,} compounds for edge generation")

    edge_buffer = []

    def ensure_node(node_id: str, name: str, kind: str):
        w_cur.execute(
            "INSERT INTO hetionet_nodes (id, name, kind) VALUES (%s,%s,%s) ON CONFLICT (id) DO NOTHING",
            (node_id, name, kind),
        )

    for row in rows:
        c_node_id = f"Compound::{row['chembl_id'] or row['name']}"
        ensure_node(c_node_id, row["name"], "Compound")

        # CtD edges
        for disease in (row["diseases"] or "").split("||"):
            disease = disease.strip()
            if disease and len(disease) > 3:
                d_node_id = f"Disease::{disease}"
                ensure_node(d_node_id, disease, "Disease")
                edge_buffer.append((c_node_id, d_node_id, "CtD"))

        # CbG edges
        for gene in (row["genes"] or "").split("||"):
            gene = gene.strip()
            if gene and len(gene) > 1:
                g_node_id = f"Gene::{gene}"
                ensure_node(g_node_id, gene, "Gene")
                edge_buffer.append((c_node_id, g_node_id, "CbG"))

    # Batch insert edges
    conn.commit()
    total = 0
    for i in range(0, len(edge_buffer), 2000):
        batch = edge_buffer[i:i+2000]
        try:
            psycopg2.extras.execute_values(
                w_cur,
                "INSERT INTO hetionet_edges (source_id, target_id, metaedge) VALUES %s ON CONFLICT DO NOTHING",
                batch,
            )
            conn.commit()
            total += len(batch)
        except Exception as e:
            conn.rollback()
            logger.warning(f"Edge batch failed: {e}")

    logger.info(f"  Generated {total:,} knowledge graph edges")
    conn.close()


def run_full_import():
    """Run the complete import pipeline."""
    from .schema import initialize_schema

    logger.info("=== NeuroRepurpose Full Import ===")
    initialize_schema()

    try:
        import_chembl_neuro_compounds()
    except psycopg2.errors.UndefinedTable as e:
        logger.warning(f"ChEMBL relational schema not present in {CHEMBL_CONN['dbname']}: {e.diag.message_primary}")
        import_chembl_snapshot_files()
    import_hetionet_data()
    import_json_drugs()
    import_json_associations_and_interactions()
    generate_chembl_edges()

    # Summary
    conn = psycopg2.connect(**NEURO_CONN)
    cur = conn.cursor()
    for table in ["compounds", "targets", "mechanisms", "indications", "hetionet_nodes", "hetionet_edges"]:
        cur.execute(f"SELECT COUNT(*) FROM {table}")
        logger.info(f"  {table}: {cur.fetchone()[0]:,} rows")
    conn.close()
    logger.info("=== Import complete ===")


if __name__ == "__main__":
    sys.path.insert(0, str(REPO_ROOT))
    run_full_import()
