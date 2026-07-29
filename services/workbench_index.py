"""
Workbench index — precompute the data-driven Repurposing Workbench landing table
═══════════════════════════════════════════════════════════════════════════════
Replaces the hand-typed RESEARCH_AREAS array on templates/index.html with a fully
computed, cached ranked table. Every column is derived from a real platform source;
nothing is hand-written. If a value cannot be computed it is left honestly empty.

Per row (all REAL, no fabrication):
  disease        MONDO preferred label
  mondo_id       MONDO id
  category       ICD10-chapter therapeutic area from the disease-value reference DB
  mechanism_genes top HGNC symbols from the SAME Open Targets disease->gene
                 resolution the scorer uses (services.disease_ontology.resolve_disease),
                 ranked by evidence-quality score. [] if the disease resolves to no genes.
  value_score    Repurposing Value Score (burden x unmet x market), precomputed in
                 services/disease_value/build_disease_table.py
  value_tier     the pharma-committee reading of that score (services.disease_value._tier)
  active_trials  live count of ACTIVE studies from ClinicalTrials.gov v2, or None if the
                 count cannot be fetched (never invented)

DISEASE LIST = the top N diseases of the real MONDO universe ranked by value_score
(ties broken deterministically by burden then label), N configurable.

Build / refresh:
    C:/Users/Soumi/AppData/Local/Python/pythoncore-3.14-64/python.exe -m services.workbench_index [N]
"""
from __future__ import annotations

import json
import logging
import sqlite3
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

_ROOT = Path(__file__).parent.parent
_VALUE_DB = _ROOT / "data" / "disease_reference" / "disease_value.db"
_OUT = _ROOT / "data" / "workbench_index.json"

DEFAULT_N = 60
# Candidate pool cap. Hundreds of orphan diseases tie at the top Repurposing Value
# score, so ranking by value alone is effectively alphabetical. We instead take a
# bounded pool of the highest-value diseases and rank THEM by live clinical-pipeline
# volume (active trials). The cap keeps the CT.gov count-fetch bounded and fast.
POOL_CAP = 300

# ClinicalTrials.gov v2 — statuses that count as an ACTIVE trial (open or ongoing).
_ACTIVE_STATUS = "RECRUITING|ACTIVE_NOT_RECRUITING|ENROLLING_BY_INVITATION|NOT_YET_RECRUITING"
_CT_STUDIES = "https://clinicaltrials.gov/api/v2/studies"


# ── Top-N disease universe (from the precomputed value DB) ───────────────────

def top_value_diseases(n: int = DEFAULT_N) -> List[Dict]:
    """Top-N diseases by Repurposing Value Score, read straight from the reference
    DB. Ties (many orphan diseases share the max score) are broken deterministically
    by burden DESC then label ASC so the table is stable across rebuilds. Each row
    carries its burden so downstream ranking can use it as a tiebreak."""
    if not _VALUE_DB.exists():
        logger.warning("workbench_index: disease_value.db missing at %s", _VALUE_DB)
        return []
    conn = sqlite3.connect(f"file:{_VALUE_DB}?mode=ro", uri=True)
    try:
        cur = conn.execute(
            "SELECT mondo_id, label, category, value_score, burden "
            "FROM diseases "
            "WHERE is_disease=1 AND value_score IS NOT NULL AND value_score > 0 "
            "ORDER BY value_score DESC, burden DESC, label ASC "
            "LIMIT ?", (int(n),))
        return [{"mondo_id": r[0], "label": r[1], "category": (r[2] or "").strip(),
                 "value_score": round(float(r[3]), 4),
                 "burden": round(float(r[4]), 4) if r[4] is not None else 0.0}
                for r in cur.fetchall()]
    finally:
        conn.close()


def value_universe_size() -> int:
    """Real count of scored diseases in the value universe (for the headline stat)."""
    if not _VALUE_DB.exists():
        return 0
    conn = sqlite3.connect(f"file:{_VALUE_DB}?mode=ro", uri=True)
    try:
        cur = conn.execute("SELECT COUNT(*) FROM diseases "
                           "WHERE is_disease=1 AND value_score IS NOT NULL AND value_score > 0")
        return int(cur.fetchone()[0])
    finally:
        conn.close()


# ── Per-disease enrichment (real sources, fail-soft to honest empty) ─────────

def mechanism_genes(disease_name: str, top_k: int = 4) -> List[str]:
    """Top HGNC symbols for the disease from the SAME Open Targets resolution the
    scorer uses, ranked by evidence-quality score. [] if it resolves to no genes."""
    try:
        from services.disease_ontology import resolve_disease
        info = resolve_disease(disease_name) or {}
    except Exception as e:
        logger.debug("workbench_index gene resolve failed for %s: %s", disease_name, e)
        return []
    targets = info.get("targets", []) or []
    ranked = sorted(
        targets,
        key=lambda t: (t.get("quality_score") or t.get("genetic_score") or t.get("score") or 0.0),
        reverse=True)
    genes: List[str] = []
    for t in ranked:
        g = (t.get("gene_symbol") or "").strip()
        if g and g not in genes:
            genes.append(g)
        if len(genes) >= top_k:
            break
    return genes


def active_trial_count(disease_name: str) -> Optional[int]:
    """Live count of ACTIVE ClinicalTrials.gov studies for the disease, or None if it
    cannot be fetched (never invented)."""
    try:
        from services import http_client
        data = http_client.get_json(
            _CT_STUDIES, default=None,
            params={"query.cond": disease_name, "filter.overallStatus": _ACTIVE_STATUS,
                    "countTotal": "true", "pageSize": 1, "format": "json"},
            timeout=25)
        if not data:
            return None
        tc = data.get("totalCount")
        return int(tc) if tc is not None else None
    except Exception as e:
        logger.debug("workbench_index trial count failed for %s: %s", disease_name, e)
        return None


# ── Build + cache ────────────────────────────────────────────────────────────

def build(n: int = DEFAULT_N, *, pool_cap: int = POOL_CAP, verbose: bool = True) -> Dict:
    """Compute the ranked table and write it to data/workbench_index.json. Re-runnable.

    Ordering (locked): the value score ties hundreds of orphan diseases at the ceiling,
    so ranking by value alone is arbitrary/alphabetical. Instead we (1) take a bounded
    POOL of the highest repurposing-value diseases, (2) fetch each one's live ACTIVE
    clinical-trial count, (3) rank the pool by active-trial volume (opportunity), and
    (4) take the top N. Mechanism genes (the expensive Open Targets call) are resolved
    ONLY for the final N, not the whole pool."""
    from services.disease_value import _tier   # the pharma-committee tier reading

    pool = top_value_diseases(pool_cap)
    if verbose:
        print(f"Pool: {len(pool)} highest-value diseases; fetching active-trial "
              f"counts (CT.gov) to rank by pipeline volume…", flush=True)

    # 1) active-trial count for every pool member (the ranking key)
    ranked: List[Dict] = []
    dropped: List[str] = []
    for i, d in enumerate(pool, 1):
        tc = active_trial_count(d["label"])
        if tc is None:
            # never rank a genuine fetch failure as 0 — drop it and note it
            dropped.append(d["label"])
        else:
            d["active_trials"] = tc
            ranked.append(d)
        if verbose and (i % 25 == 0 or i == len(pool)):
            print(f"  trial counts {i}/{len(pool)} (kept {len(ranked)}, dropped {len(dropped)})",
                  flush=True)

    # 2) rank by active trials DESC, deterministic tiebreak value/burden/label
    ranked.sort(key=lambda d: (-d["active_trials"], -d["value_score"],
                               -d.get("burden", 0.0), d["label"]))
    final = ranked[:int(n)]

    # 3) resolve mechanism genes only for the final N
    rows: List[Dict] = []
    for i, d in enumerate(final, 1):
        name = d["label"]
        genes = mechanism_genes(name)
        rows.append({
            "disease": name,
            "mondo_id": d["mondo_id"],
            "category": d["category"] or "Unclassified",
            "mechanism_genes": genes,
            "value_score": d["value_score"],
            "value_tier": _tier(d["value_score"]),
            "active_trials": d["active_trials"],
        })
        if verbose:
            print(f"[{i}/{len(final)}] {name[:38]:<38} trials={d['active_trials']:<4} "
                  f"tier={_tier(d['value_score'])} genes={genes}", flush=True)

    payload = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_requested": int(n),
        "n_rows": len(rows),
        "pool_size": len(pool),
        "pool_cap": int(pool_cap),
        "dropped_no_trial_count": dropped,
        "universe_size": value_universe_size(),
        "ranking": "within top repurposing-value tier, ordered by active clinical trial volume",
        "sources": {
            "value_score": "services/disease_value/build_disease_table.py (MONDO burden x unmet x market)",
            "mechanism_genes": "services/disease_ontology.py resolve_disease (Open Targets associatedTargets)",
            "category": "disease_value.db category (ICD10-chapter therapeutic area)",
            "active_trials": "ClinicalTrials.gov API v2 studies countTotal, active statuses",
        },
        "rows": rows,
    }
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    _OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    if verbose:
        with_genes = sum(1 for r in rows if r["mechanism_genes"])
        print(f"\nDONE -> {_OUT}\n  rows={len(rows)}  with_genes={with_genes}  "
              f"pool={len(pool)}  dropped_no_count={len(dropped)}  "
              f"universe={payload['universe_size']}", flush=True)
    return payload


def load() -> Optional[Dict]:
    """Return the cached workbench index, or None if it has not been built."""
    if not _OUT.exists():
        return None
    try:
        return json.loads(_OUT.read_text(encoding="utf-8"))
    except Exception as e:
        logger.warning("workbench_index: cache unreadable: %s", e)
        return None


def fallback_rows(n: int = DEFAULT_N) -> Dict:
    """Honest fallback when the cache is missing: the real top-N value diseases with
    category and value tier from the DB, but NO invented mechanisms or trial counts
    (mechanism blank, trials null). Never returns the old hardcoded array."""
    from services.disease_value import _tier
    diseases = top_value_diseases(n)
    rows = [{
        "disease": d["label"], "mondo_id": d["mondo_id"],
        "category": d["category"] or "Unclassified",
        "mechanism_genes": [], "value_score": d["value_score"],
        "value_tier": _tier(d["value_score"]), "active_trials": None,
    } for d in diseases]
    return {"built_at": None, "n_requested": int(n), "n_rows": len(rows),
            "universe_size": value_universe_size(), "computing": True, "rows": rows}


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_N
    t0 = time.time()
    build(n)
    print(f"elapsed {time.time() - t0:.1f}s", flush=True)
