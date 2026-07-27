"""
Mechanism-Direction Concordance  -  ChEMBL 33  vs  IUPHAR/GtoPdb
═══════════════════════════════════════════════════════════════════════════
Independent-source accuracy check for DRUG MECHANISM DIRECTION — the input the
platform's direction-aware pathway scoring depends on.

The platform reads each drug's mechanism action_type from ChEMBL (INHIBITOR /
AGONIST / ANTAGONIST / …) and converts it to a signed perturbation direction
(inhibit = −1, activate = +1). If that direction is wrong, the direction-aware
pathway score (and the "suppresses vs activates pathway" verdict on the Pathway
Screen) is wrong. So we check it against a second, independently-curated source.

Method:
  1. GtoPdb interactions give (ligand, target gene symbol, interaction Type).
  2. Map GtoPdb ligand -> ChEMBL ID (from GtoPdb ligands.csv).
  3. ChEMBL drug_mechanism gives (molecule, target gene symbol, action_type).
  4. For (molecule, gene) pairs present in BOTH, map each side to a sign
     (−1 inhibit/antagonise/block, +1 activate/agonise) and check agreement.

Read-only against chembl_33. GtoPdb files default to data/external/gtopdb/.
Writes validation/mechanisms_results.json + run log.

Usage:
    python validation/validate_mechanisms.py [gtp_interactions.csv] [gtp_ligands.csv]
"""

import sys
import csv
import json
import datetime
from pathlib import Path

from sqlalchemy import text

try:                                   # Windows consoles default to cp1252
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from validation.validate_concordance import get_engine, _csv_rows  # reuse infra

DEFAULT_INTER = ROOT / "data" / "external" / "gtopdb" / "interactions.csv"
DEFAULT_LIGS  = ROOT / "data" / "external" / "gtopdb" / "ligands.csv"

# Direction sign mapping (substring match -> robust to GtoPdb / ChEMBL wording).
_NEG = ("INHIBIT", "ANTAGONIST", "ANTAGONISM", "BLOCK", "NEGATIVE", "INVERSE AGONIST",
        "SUPPRESS", "DEGRAD", "DISRUPT", "ANTISENSE")
_POS = ("AGONIST", "ACTIVAT", "POSITIVE", "OPENER", "INDUC", "STABILIS", "STABILIZ")


def _sign(label: str):
    """Map an interaction/action label to −1, +1, or None (no asserted direction)."""
    s = (label or "").strip().upper()
    if not s:
        return None
    # inverse agonist / negative modulator must beat the bare 'AGONIST'/'POSITIVE' check
    if "INVERSE AGONIST" in s or "NEGATIVE" in s:
        return -1
    if any(t in s for t in _NEG):
        return -1
    if any(t in s for t in _POS):
        return 1
    return None


def run(interactions_path, ligands_path):
    log_lines = []

    def log(msg=""):
        print(msg)
        log_lines.append(str(msg))

    log("=" * 72)
    log("MECHANISM-DIRECTION CONCORDANCE  -  ChEMBL 33  vs  IUPHAR/GtoPdb")
    log(f"Run: {datetime.datetime.now().isoformat(timespec='seconds')}")
    log("=" * 72)

    # 1. GtoPdb ligand_id -> ChEMBL ID
    lig_to_chembl = {}
    for row in _csv_rows(ligands_path):
        lid = (row.get("Ligand ID") or "").strip()
        ch = (row.get("ChEMBL ID") or "").strip()
        if lid and ch:
            lig_to_chembl[lid] = ch
    log(f"\nGtoPdb ligands with a ChEMBL ID: {len(lig_to_chembl):,}")

    # 2. GtoPdb directional interactions -> {(chembl_id, GENE): sign}
    gtp_dir = {}
    gtp_rows = 0
    for row in _csv_rows(interactions_path):
        if (row.get("Target Species") or "").strip() != "Human":
            continue
        gene = (row.get("Target Gene Symbol") or "").strip().upper()
        lid = (row.get("Ligand ID") or "").strip()
        ch = lig_to_chembl.get(lid)
        if not (gene and ch):
            continue
        sign = _sign(row.get("Type"))
        if sign is None:
            sign = _sign(row.get("Action"))
        if sign is None:
            continue
        gtp_rows += 1
        key = (ch, gene)
        # keep the first asserted direction; conflicts within GtoPdb are rare
        gtp_dir.setdefault(key, sign)
    log(f"GtoPdb human directional (ligand,gene) interactions: {gtp_rows:,} "
        f"-> {len(gtp_dir):,} unique (drug,gene) pairs")

    chembl_ids = sorted({k[0] for k in gtp_dir})

    # 3. ChEMBL drug_mechanism direction -> {(chembl_id, GENE): sign}
    eng = get_engine()
    chembl_dir = {}
    with eng.connect() as c:
        rows = c.execute(text(
            """
            SELECT md.chembl_id, UPPER(csyn.component_synonym), dm.action_type
            FROM molecule_dictionary md
            JOIN drug_mechanism dm        ON dm.molregno = md.molregno
            JOIN target_components tc     ON tc.tid = dm.tid
            JOIN component_synonyms csyn  ON csyn.component_id = tc.component_id
                                         AND csyn.syn_type = 'GENE_SYMBOL'
            WHERE md.chembl_id = ANY(:ids) AND dm.action_type IS NOT NULL
            """), {"ids": chembl_ids}).fetchall()
    for ch, gene, action in rows:
        sign = _sign(action)
        if sign is not None:
            chembl_dir.setdefault((ch, gene), sign)
    log(f"ChEMBL directional (drug,gene) mechanism pairs (for these drugs): {len(chembl_dir):,}")

    # 4. Compare on the intersection
    shared = sorted(set(gtp_dir) & set(chembl_dir))
    agree = 0
    disagreements = []
    for key in shared:
        if gtp_dir[key] == chembl_dir[key]:
            agree += 1
        elif len(disagreements) < 25:
            disagreements.append({"chembl_id": key[0], "gene": key[1],
                                  "chembl_sign": chembl_dir[key], "gtopdb_sign": gtp_dir[key]})
    compared = len(shared)
    pct = round(100.0 * agree / compared, 1) if compared else 0.0

    log(f"\n(drug,gene) pairs directional in BOTH sources: {compared:,}")
    log(f"Direction AGREEMENT: {agree:,}/{compared:,}  =  {pct}%")
    if disagreements:
        log("Sample disagreements (sign convention: -1 inhibit/antagonise, +1 activate):")
        for d in disagreements[:8]:
            log(f"  {d['chembl_id']:<16} {d['gene']:<10} ChEMBL={d['chembl_sign']:+d}  GtoPdb={d['gtopdb_sign']:+d}")

    severity = "INFO" if pct >= 90 else ("WARN" if pct >= 80 else "FAIL")
    result = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "comparison": "ChEMBL 33 mechanism action_type vs IUPHAR/GtoPdb interaction Type",
        "validates": "Direction sign used by the direction-aware pathway score and the "
                     "Pathway Screen 'suppresses vs activates' verdict.",
        "pairs_compared": compared,
        "direction_agreement_pct": pct,
        "agree": agree,
        "gtopdb_pairs": len(gtp_dir),
        "chembl_pairs": len(chembl_dir),
        "disagreement_examples": disagreements,
        "findings": [{
            "id": "MECH-01",
            "severity": severity,
            "title": "Mechanism-direction concordance with IUPHAR/GtoPdb",
            "detail": f"Of {compared:,} (drug, target-gene) pairs whose interaction direction is "
                      f"asserted in BOTH ChEMBL and the independently-curated IUPHAR/GtoPdb, "
                      f"{pct}% assign the SAME direction (inhibit vs activate). This validates the "
                      f"signed mechanism input the direction-aware pathway score relies on.",
            "capa": "Pairs with opposite signs are flagged for review; the platform keeps ChEMBL's "
                    "curated action_type but the disagreement list is the audit trail.",
        }],
    }
    (HERE / "mechanisms_results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (HERE / "mechanisms_run.log").write_text("\n".join(log_lines), encoding="utf-8")
    log(f"\n  Wrote: {HERE / 'mechanisms_results.json'}")
    return result


if __name__ == "__main__":
    inter = sys.argv[1] if len(sys.argv) > 1 else str(DEFAULT_INTER)
    ligs = sys.argv[2] if len(sys.argv) > 2 else str(DEFAULT_LIGS)
    try:
        run(inter, ligs)
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(1)
