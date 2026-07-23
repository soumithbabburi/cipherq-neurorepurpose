"""
MeSH Disease Ontology Importer
Fetches all 1,938 unique MeSH IDs from ChEMBL drug_indication, resolves full
descriptor data (tree numbers, entry terms, parent IDs, scope notes) via the
NLM MeSH SPARQL endpoint (using curl since Python 3.14 TLS rejects NLM's server),
and stores in the mesh_diseases table.

Also:
- Backfills indications.mesh_id from ChEMBL
- Imports pChEMBL activity data for relevant target/compound pairs
"""

import json
import logging
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import quote

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


CHEMBL_CONN = _conn_config("CHEMBL_DB", "chembl_33")
NEURO_CONN  = {
    "host": os.getenv("DB_HOST", CHEMBL_CONN["host"]),
    "port": int(os.getenv("DB_PORT", str(CHEMBL_CONN["port"]))),
    "user": os.getenv("DB_USER", CHEMBL_CONN["user"]),
    "password": os.getenv("DB_PASSWORD", CHEMBL_CONN["password"]),
    "dbname": os.getenv("DB_NAME", "neurorepurpose"),
}

SPARQL_ENDPOINT = "https://id.nlm.nih.gov/mesh/sparql"
MESH_BASE       = "http://id.nlm.nih.gov/mesh/"
VOCAB           = "http://id.nlm.nih.gov/mesh/vocab#"


# ─── curl-based HTTP helper ───────────────────────────────────────────────────

def _curl_get(url: str, params: dict = None, retries: int = 3) -> Optional[dict]:
    """Fetch a URL via curl (bypasses Python 3.14 TLS issues with NLM)."""
    if params:
        query_str = "&".join(f"{k}={quote(str(v))}" for k, v in params.items())
        url = f"{url}?{query_str}"

    for attempt in range(retries):
        try:
            result = subprocess.run(
                ["curl", "-s", "-L", "--max-time", "30",
                 "-H", "Accept: application/sparql-results+json, application/json",
                 "-H", "User-Agent: NeuroRepurpose/1.0 (research)",
                 url],
                capture_output=True, text=True, timeout=35,
            )
            if result.returncode == 0 and result.stdout.strip():
                return json.loads(result.stdout)
        except (subprocess.TimeoutExpired, json.JSONDecodeError) as e:
            logger.warning(f"curl attempt {attempt+1} failed for {url[:80]}: {e}")
            time.sleep(1.5 * (attempt + 1))
    return None


def _sparql_query(query: str) -> List[dict]:
    """Run a SPARQL query against the NLM endpoint, return bindings list."""
    data = _curl_get(SPARQL_ENDPOINT, {"query": query, "format": "JSON", "limit": "500"})
    if data is None:
        return []
    return data.get("results", {}).get("bindings", [])


def _val(binding: dict, key: str) -> str:
    return binding.get(key, {}).get("value", "")


# ─── Fetch all MeSH IDs from ChEMBL ─────────────────────────────────────────

def _get_chembl_mesh_ids() -> List[tuple]:
    """Return all (mesh_id, mesh_heading) pairs from ChEMBL drug_indication."""
    conn = psycopg2.connect(**CHEMBL_CONN)
    cur  = conn.cursor()
    cur.execute(
        "SELECT DISTINCT mesh_id, mesh_heading "
        "FROM drug_indication "
        "WHERE mesh_id IS NOT NULL AND mesh_id != '' "
        "ORDER BY mesh_id"
    )
    rows = cur.fetchall()
    conn.close()
    logger.info(f"Found {len(rows):,} unique MeSH IDs in ChEMBL")
    return rows


# ─── Fetch descriptor data in batches via SPARQL ─────────────────────────────

def _batch_fetch_descriptors(mesh_ids: List[str], batch_size: int = 50) -> Dict[str, dict]:
    """
    Fetch tree numbers and broader (parent) descriptors for a list of MeSH IDs.
    Returns dict: mesh_id -> {tree_numbers: [...], parent_ids: [...]}
    """
    results = {}
    total   = len(mesh_ids)

    for i in range(0, total, batch_size):
        batch = mesh_ids[i : i + batch_size]
        values = " ".join(f"<{MESH_BASE}{mid}>" for mid in batch)

        query = f"""
        PREFIX meshv: <{VOCAB}>
        PREFIX rdfs:  <http://www.w3.org/2000/01/rdf-schema#>

        SELECT ?d ?tree ?broader WHERE {{
          VALUES ?d {{ {values} }}
          OPTIONAL {{ ?d meshv:treeNumber ?tree }}
          OPTIONAL {{ ?d meshv:broaderDescriptor ?broader }}
        }}
        """

        bindings = _sparql_query(query)
        for b in bindings:
            mid = _val(b, "d").replace(MESH_BASE, "").split("/")[-1]
            if mid not in results:
                results[mid] = {"tree_numbers": [], "parent_ids": []}
            tree = _val(b, "tree")
            if tree:
                t = tree.replace(MESH_BASE, "").replace(MESH_BASE.rstrip("/") + "/", "")
                # tree numbers look like http://id.nlm.nih.gov/mesh/C10.228... — extract suffix
                for seg in [MESH_BASE + "2025/", MESH_BASE + "2024/", MESH_BASE + "2023/", MESH_BASE]:
                    t = t.replace(seg, "")
                if t and t not in results[mid]["tree_numbers"]:
                    results[mid]["tree_numbers"].append(t)
            broader = _val(b, "broader")
            if broader:
                pid = broader.split("/")[-1]
                if pid and pid not in results[mid]["parent_ids"]:
                    results[mid]["parent_ids"].append(pid)

        done = min(i + batch_size, total)
        logger.info(f"  Descriptor batch {done}/{total} — resolved {len(results):,}")
        time.sleep(0.3)

    return results


def _batch_fetch_entry_terms(mesh_ids: List[str], batch_size: int = 30) -> Dict[str, List[str]]:
    """
    Fetch entry terms (synonyms/aliases) for each MeSH ID via SPARQL.
    Returns dict: mesh_id -> [term1, term2, ...]
    """
    results = {}
    total   = len(mesh_ids)

    for i in range(0, total, batch_size):
        batch  = mesh_ids[i : i + batch_size]
        values = " ".join(f"<{MESH_BASE}{mid}>" for mid in batch)

        # MeSH descriptors expose their synonyms via (preferred|narrower) Concept →
        # Term → (prefLabel|altLabel). The earlier `meshv:concept`/`rdfs:label` path
        # returns nothing on the live NLM endpoint, which left entry_terms empty.
        query = f"""
        PREFIX meshv: <{VOCAB}>
        PREFIX rdfs:  <http://www.w3.org/2000/01/rdf-schema#>

        SELECT ?d ?term WHERE {{
          VALUES ?d {{ {values} }}
          ?d meshv:concept|meshv:preferredConcept ?concept .
          ?concept meshv:term ?t .
          ?t meshv:prefLabel|meshv:altLabel ?term .
        }}
        """

        bindings = _sparql_query(query)
        for b in bindings:
            mid  = _val(b, "d").split("/")[-1]
            term = _val(b, "term")
            if mid and term:
                results.setdefault(mid, [])
                if term not in results[mid]:
                    results[mid].append(term)

        done = min(i + batch_size, total)
        if done % 150 == 0 or done == total:
            logger.info(f"  Entry-term batch {done}/{total}")
        time.sleep(0.3)

    return results


# ─── Import into DB ───────────────────────────────────────────────────────────

def import_mesh_diseases():
    """Main MeSH import: fetch from NLM SPARQL + store in mesh_diseases table."""
    logger.info("=== MeSH Disease Ontology Import ===")

    rows       = _get_chembl_mesh_ids()
    mesh_ids   = [r[0] for r in rows]
    heading_map = {r[0]: r[1] for r in rows}

    logger.info("Fetching tree numbers and parent IDs...")
    descriptor_data = _batch_fetch_descriptors(mesh_ids)

    logger.info("Fetching entry terms (synonyms)...")
    entry_term_data = _batch_fetch_entry_terms(mesh_ids)

    conn = psycopg2.connect(**NEURO_CONN)
    cur  = conn.cursor()
    inserted = updated = 0

    for mid in mesh_ids:
        heading     = heading_map.get(mid, mid)
        desc        = descriptor_data.get(mid, {})
        tree_nums   = desc.get("tree_numbers", [])
        parent_ids  = desc.get("parent_ids", [])
        entry_terms = entry_term_data.get(mid, [])
        # Always include the heading itself as an entry term
        if heading and heading not in entry_terms:
            entry_terms.insert(0, heading)

        cur.execute(
            """
            INSERT INTO mesh_diseases (mesh_id, heading, tree_numbers, entry_terms, parent_ids)
            VALUES (%s, %s, %s, %s, %s)
            ON CONFLICT (mesh_id) DO UPDATE
                SET heading     = EXCLUDED.heading,
                    tree_numbers = EXCLUDED.tree_numbers,
                    entry_terms  = EXCLUDED.entry_terms,
                    parent_ids   = EXCLUDED.parent_ids
            """,
            (mid, heading, tree_nums, entry_terms, parent_ids),
        )
        inserted += 1

    conn.commit()
    logger.info(f"Upserted {inserted:,} MeSH disease records")

    # Backfill indications.mesh_id from ChEMBL
    logger.info("Backfilling indications.mesh_id from ChEMBL...")
    src = psycopg2.connect(**CHEMBL_CONN)
    src_cur = src.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    src_cur.execute(
        "SELECT DISTINCT mesh_id, mesh_heading FROM drug_indication "
        "WHERE mesh_id IS NOT NULL AND mesh_id != ''"
    )
    mesh_map = {r["mesh_heading"].lower(): r["mesh_id"] for r in src_cur.fetchall()}
    src.close()

    cur.execute("SELECT id, disease FROM indications WHERE mesh_id IS NULL")
    to_update = cur.fetchall()
    updated = 0
    for row_id, disease in to_update:
        mid = mesh_map.get((disease or "").lower())
        if mid:
            cur.execute("UPDATE indications SET mesh_id = %s WHERE id = %s", (mid, row_id))
            updated += 1

    conn.commit()
    conn.close()
    logger.info(f"Backfilled mesh_id on {updated:,} indication rows")
    logger.info("=== MeSH import complete ===")


# ─── Activity data import ─────────────────────────────────────────────────────

def import_activity_data():
    """
    Import pChEMBL Ki/IC50 values for our neuro compounds at their known targets.
    Only imports records with pchembl_value (quantitative, curated data).
    """
    logger.info("=== Activity Data Import ===")

    conn = psycopg2.connect(**NEURO_CONN)
    cur  = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    # Get all our compound ChEMBL IDs → internal IDs
    cur.execute("SELECT id, chembl_id FROM compounds WHERE chembl_id IS NOT NULL")
    compound_map = {r["chembl_id"]: r["id"] for r in cur.fetchall()}

    # Get all our target ChEMBL IDs → internal IDs
    cur.execute("SELECT id, chembl_tid FROM targets WHERE chembl_tid IS NOT NULL")
    target_map = {r["chembl_tid"]: r["id"] for r in cur.fetchall()}

    if not compound_map or not target_map:
        logger.warning("No compounds or targets in DB — skipping activity import")
        conn.close()
        return

    chembl_cids = list(compound_map.keys())
    chembl_tids = list(target_map.keys())

    src = psycopg2.connect(**CHEMBL_CONN)
    src_cur = src.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    cp_placeholders = ",".join(["%s"] * len(chembl_cids))
    td_placeholders = ",".join(["%s"] * len(chembl_tids))

    src_cur.execute(
        f"""
        SELECT md.chembl_id AS compound_chembl,
               td.chembl_id AS target_chembl,
               a.standard_type AS activity_type,
               a.pchembl_value,
               a.standard_value,
               a.standard_units
        FROM activities a
        JOIN assays ass ON ass.assay_id = a.assay_id
        JOIN target_dictionary td ON td.tid = ass.tid
        JOIN molecule_dictionary md ON md.molregno = a.molregno
        WHERE md.chembl_id IN ({cp_placeholders})
          AND td.chembl_id IN ({td_placeholders})
          AND a.pchembl_value IS NOT NULL
          AND a.standard_type IN ('Ki', 'IC50', 'Kd', 'EC50')
          AND (a.data_validity_comment IS NULL OR a.data_validity_comment = 'Manually validated')
        ORDER BY a.pchembl_value DESC
        """,
        chembl_cids + chembl_tids,
    )
    rows = src_cur.fetchall()
    src.close()
    logger.info(f"  Got {len(rows):,} activity records from ChEMBL")

    w_cur = conn.cursor()
    # Clear and re-insert
    w_cur.execute("DELETE FROM compound_activities")
    batch = []
    for r in rows:
        cid = compound_map.get(r["compound_chembl"])
        tid = target_map.get(r["target_chembl"])
        if cid and tid:
            batch.append((
                cid, tid,
                r["activity_type"],
                float(r["pchembl_value"]) if r["pchembl_value"] else None,
                float(r["standard_value"]) if r["standard_value"] else None,
                r["standard_units"],
            ))

    if batch:
        psycopg2.extras.execute_values(
            w_cur,
            "INSERT INTO compound_activities (compound_id, target_id, activity_type, pchembl_value, standard_value, standard_units) VALUES %s",
            batch, page_size=1000,
        )

    conn.commit()
    conn.close()
    logger.info(f"  Imported {len(batch):,} activity records")
    logger.info("=== Activity import complete ===")


if __name__ == "__main__":
    sys.path.insert(0, str(REPO_ROOT))
    import_mesh_diseases()
    import_activity_data()
