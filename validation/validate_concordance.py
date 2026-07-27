"""
Cross-Database Concordance  -  ChEMBL 33  vs  Guide to Pharmacology (IUPHAR)
═══════════════════════════════════════════════════════════════════════════
Independent-source accuracy check. For compound-target affinity pairs that
appear in BOTH ChEMBL and the IUPHAR/Guide to Pharmacology (GtoPdb) - a
separately-curated gold-standard database - do the two agree?

This is the "accuracy" dimension of the data-quality scorecard: not ChEMBL
checking itself, but ChEMBL measured against a second, independent authority.

Method:
  1. GtoPdb interactions give (ligand, target UniProt, pKi/pIC50/pEC50/pKd value).
  2. Map GtoPdb ligand -> ChEMBL ID (from GtoPdb ligands.csv).
  3. Map UniProt -> ChEMBL target (component_sequences.accession -> tid).
  4. ChEMBL median pchembl for the same molecule+target+type (high-confidence).
  5. Concordant if |ChEMBL_median - GtoPdb_value| <= 1.0 log unit.

Read-only against chembl_33. GtoPdb files are downloaded to /tmp by the caller
(or pass paths). Writes validation/concordance_results.json + run log.

Usage:
    python validation/validate_concordance.py [gtp_interactions.csv] [gtp_ligands.csv]
"""

import os
import sys
import csv
import json
import datetime
import statistics
from pathlib import Path

from sqlalchemy import create_engine, text
from sqlalchemy.engine import URL

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

CONCORDANCE_TOLERANCE_LOG = 1.0          # within 1 log = concordant
HIGH_CONFIDENCE = 8
# GtoPdb pUnit -> ChEMBL standard_type
PUNIT_TO_TYPE = {"pKi": "Ki", "pIC50": "IC50", "pEC50": "EC50", "pKd": "Kd"}


def load_env():
    env = {}
    envfile = ROOT / ".env"
    if envfile.exists():
        for line in envfile.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                env[k.strip()] = v.strip()
    return env


def get_engine():
    env = load_env()

    def pick(*keys, default=""):
        for k in keys:
            v = env.get(k) or os.environ.get(k)
            if v:
                return v
        return default

    uri = URL.create(
        "postgresql+psycopg2",
        username=pick("DB_USER", "CHEMBL_DB_USER", default="babburisoumith"),
        password=pick("DB_PASSWORD", "CHEMBL_DB_PASSWORD", "PGPASSWORD") or None,
        host=pick("DB_HOST", "CHEMBL_DB_HOST", default="localhost"),
        port=int(pick("DB_PORT", "CHEMBL_DB_PORT", default="5433")),
        database="chembl_33",
    )
    return create_engine(uri, connect_args={"connect_timeout": 30})


def _csv_rows(path, skip_comment=True):
    """Yield dict rows; GtoPdb files have a version-comment first line."""
    with open(path, encoding="utf-8") as f:
        if skip_comment:
            next(f)  # version comment
        for row in csv.DictReader(f):
            yield row


def run(interactions_path, ligands_path):
    log_lines = []

    def log(msg=""):
        print(msg)
        log_lines.append(str(msg))

    log("=" * 72)
    log("CROSS-DB CONCORDANCE  -  ChEMBL 33  vs  IUPHAR/GtoPdb")
    log(f"Run: {datetime.datetime.now().isoformat(timespec='seconds')}")
    log("=" * 72)

    # 1. ligand_id -> ChEMBL ID
    lig_to_chembl = {}
    for row in _csv_rows(ligands_path):
        lid = (row.get("Ligand ID") or "").strip()
        ch = (row.get("ChEMBL ID") or "").strip()
        if lid and ch:
            lig_to_chembl[lid] = ch
    log(f"\nGtoPdb ligands with a ChEMBL ID: {len(lig_to_chembl):,}")

    # 2. external pairs (chembl_id, uniprot, std_type, ext_pval)
    ext_pairs = []
    for row in _csv_rows(interactions_path):
        punit = (row.get("Affinity Units") or "").strip()
        std_type = PUNIT_TO_TYPE.get(punit)
        if not std_type:
            continue
        if (row.get("Original Affinity Relation") or "").strip() not in ("=", ""):
            continue
        med = (row.get("Affinity Median") or "").strip()
        uni = (row.get("Target UniProt ID") or "").strip()
        lid = (row.get("Ligand ID") or "").strip()
        ch = lig_to_chembl.get(lid)
        if not (med and uni and ch):
            continue
        # UniProt field may carry a complex; take the first accession token.
        uni = uni.replace("|", " ").split()[0]
        try:
            ext_pairs.append((ch, uni, std_type, float(med)))
        except ValueError:
            continue
    log(f"GtoPdb affinity pairs (mappable, relation '='): {len(ext_pairs):,}")

    chembl_ids = sorted({p[0] for p in ext_pairs})
    uniprots = sorted({p[1] for p in ext_pairs})

    eng = get_engine()

    def rows(sql, **p):
        with eng.connect() as c:
            return c.execute(text(sql), p).fetchall()

    # 3. ChEMBL maps
    mol_map = {r[0]: r[1] for r in rows(
        "SELECT chembl_id, molregno FROM molecule_dictionary WHERE chembl_id = ANY(:ids)",
        ids=chembl_ids)}
    uni_map = {}
    for r in rows(
        "SELECT cs.accession, tc.tid FROM component_sequences cs "
        "JOIN target_components tc ON tc.component_id = cs.component_id "
        "WHERE cs.accession = ANY(:accs)", accs=uniprots):
        uni_map.setdefault(r[0], set()).add(r[1])
    log(f"ChEMBL molecule matches: {len(mol_map):,} / {len(chembl_ids):,} ligands")
    log(f"ChEMBL target matches:   {len(uni_map):,} / {len(uniprots):,} UniProt IDs")

    # 4. ChEMBL median pchembl per (molregno, tid, type) - high-confidence
    molregnos = sorted(set(mol_map.values()))
    chembl_med = {}
    if molregnos:
        for r in rows(
            """
            SELECT a.molregno, ass.tid, a.standard_type,
                   percentile_cont(0.5) WITHIN GROUP (ORDER BY a.pchembl_value) AS med
            FROM activities a JOIN assays ass ON a.assay_id = ass.assay_id
            WHERE a.molregno = ANY(:mols)
              AND a.standard_type = ANY(:types)
              AND a.pchembl_value IS NOT NULL
              AND ass.confidence_score >= :cs
              AND a.data_validity_comment IS NULL
              AND a.standard_relation = '='
            GROUP BY a.molregno, ass.tid, a.standard_type
            """,
            mols=molregnos, types=list(PUNIT_TO_TYPE.values()), cs=HIGH_CONFIDENCE):
            chembl_med[(r[0], r[1], r[2])] = float(r[3])
    log(f"ChEMBL high-confidence median values indexed: {len(chembl_med):,}")

    # 5. compare
    per_type = {t: {"compared": 0, "concordant": 0, "absdiff": []} for t in PUNIT_TO_TYPE.values()}
    all_absdiff = []
    examples = []
    for ch, uni, std_type, ext_pval in ext_pairs:
        molregno = mol_map.get(ch)
        tids = uni_map.get(uni)
        if molregno is None or not tids:
            continue
        # best (closest) matching ChEMBL median across the target's tids
        best = None
        for tid in tids:
            cm = chembl_med.get((molregno, tid, std_type))
            if cm is not None and (best is None or abs(cm - ext_pval) < abs(best - ext_pval)):
                best = cm
        if best is None:
            continue
        diff = abs(best - ext_pval)
        per_type[std_type]["compared"] += 1
        per_type[std_type]["absdiff"].append(diff)
        all_absdiff.append(diff)
        if diff <= CONCORDANCE_TOLERANCE_LOG:
            per_type[std_type]["concordant"] += 1
        if len(examples) < 20:
            examples.append({"chembl_id": ch, "uniprot": uni, "type": std_type,
                             "chembl": round(best, 2), "gtopdb": round(ext_pval, 2),
                             "abs_diff": round(diff, 2)})

    # aggregate
    by_type = {}
    for t, d in per_type.items():
        if d["compared"]:
            by_type[t] = {
                "compared": d["compared"],
                "concordant": d["concordant"],
                "concordance_pct": round(100.0 * d["concordant"] / d["compared"], 1),
                "median_abs_diff_log": round(statistics.median(d["absdiff"]), 3),
            }
    total_compared = sum(d["compared"] for d in per_type.values())
    total_conc = sum(d["concordant"] for d in per_type.values())
    overall_pct = round(100.0 * total_conc / total_compared, 1) if total_compared else 0.0
    median_abs = round(statistics.median(all_absdiff), 3) if all_absdiff else None

    log("\nCONCORDANCE (within 1.0 log of GtoPdb):")
    for t, d in by_type.items():
        log(f"  {t:5s}: {d['concordant']:>5,}/{d['compared']:<5,}  {d['concordance_pct']:>5}%  "
            f"median abs diff={d['median_abs_diff_log']} log")
    log(f"  ----")
    log(f"  ALL  : {total_conc:,}/{total_compared:,}  {overall_pct}%  median abs diff={median_abs} log")

    result = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "comparison": "ChEMBL 33 vs IUPHAR/Guide to Pharmacology",
        "tolerance_log": CONCORDANCE_TOLERANCE_LOG,
        "pairs_compared": total_compared,
        "overall_concordance_pct": overall_pct,
        "overall_median_abs_diff_log": median_abs,
        "by_type": by_type,
        "examples": examples,
        "findings": [{
            "id": "CONC-01",
            "severity": "INFO",
            "title": "Cross-database concordance with IUPHAR/GtoPdb",
            "detail": f"Of {total_compared:,} compound-target affinity pairs present in BOTH ChEMBL and "
                      f"the independently-curated IUPHAR/GtoPdb, {overall_pct}% agree within "
                      f"{CONCORDANCE_TOLERANCE_LOG:.0f} log unit (median absolute difference {median_abs} log). "
                      f"This is an external accuracy check - ChEMBL measured against a second authority, not itself.",
            "capa": "Pairs disagreeing by >1 log are candidates for review/exclusion from the certified tier.",
        }],
    }
    (HERE / "concordance_results.json").write_text(json.dumps(result, indent=2))
    (HERE / "concordance_run.log").write_text("\n".join(log_lines))
    log(f"\n  Wrote: {HERE / 'concordance_results.json'}")
    return result


if __name__ == "__main__":
    inter = sys.argv[1] if len(sys.argv) > 1 else "/tmp/gtp_interactions.csv"
    ligs = sys.argv[2] if len(sys.argv) > 2 else "/tmp/gtp_ligands.csv"
    try:
        run(inter, ligs)
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(1)
