"""
Disease Entity-Resolution Audit  -  MeSH mapping correctness
════════════════════════════════════════════════════════════
The highest-risk silent error in the platform: a disease name resolving to the
WRONG MeSH descriptor. Every downstream candidate would then be garbage while
looking fine. This audit measures resolution accuracy by round-trip identity:

  - Heading round-trip: a descriptor's own heading must resolve back to itself.
  - Synonym round-trip: a descriptor's entry term (synonym) must resolve to the
    SAME descriptor (the real test - synonyms are where mis-mapping happens).
  - Alias resolution: common clinical abbreviations resolve to the right disease.

Read-only. Writes validation/resolution_results.json + run log.

Usage:
    python validation/validate_resolution.py
"""

import os
import sys
import json
import random
import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import psycopg2
import psycopg2.extras
from config import db_params
from services.disease_resolver import resolve_disease

SAMPLE_N = 300
random.seed(17)

# Hand-curated GOLDEN SET - tricky disease synonyms / abbreviations / multi-word
# clinical names paired with their KNOWN-correct MeSH descriptor (expected heading
# substring). Unlike round-trip identity (which proves CODE CONSISTENCY), this
# tests CLINICAL-MAPPING CORRECTNESS against externally-known ground truth.
# Scoped to diseases that are in the loaded MeSH subset so it measures correctness,
# not coverage.
GOLDEN_SET = {
    # abbreviations
    "PD": "parkinson", "AD": "alzheimer", "ALS": "amyotrophic lateral sclerosis",
    "MS": "multiple sclerosis", "HD": "huntington", "FTD": "frontotemporal",
    "ADHD": "attention deficit", "PTSD": "stress disorders, post-traumatic",
    "OCD": "obsessive-compulsive", "TBI": "brain injuries",
    "T2D": "diabetes mellitus, type 2", "T2DM": "diabetes mellitus, type 2",
    "RA": "arthritis, rheumatoid", "SLE": "lupus erythematosus, systemic",
    "COPD": "pulmonary disease, chronic obstructive", "IPF": "idiopathic pulmonary fibrosis",
    "NSCLC": "carcinoma, non-small-cell lung", "AML": "leukemia, myeloid, acute",
    "HCC": "carcinoma, hepatocellular", "GBM": "glioblastoma", "MM": "multiple myeloma",
    "CKD": "renal insufficiency, chronic", "PAH": "hypertension, pulmonary",
    "AFib": "atrial fibrillation", "CAD": "coronary artery disease",
    "SMA": "muscular atrophy, spinal", "DMD": "muscular dystrophy, duchenne",
    "NMO": "neuromyelitis optica", "NASH": "non-alcoholic fatty liver",
    # multi-word names / lay synonyms
    "Lou Gehrig's disease": "amyotrophic lateral sclerosis",
    "rheumatoid arthritis": "arthritis, rheumatoid",
    "non-small cell lung cancer": "carcinoma, non-small-cell lung",
    "acute myeloid leukemia": "leukemia, myeloid, acute",
    "type 2 diabetes": "diabetes mellitus, type 2",
}


def _split_entry_terms(val):
    if val is None:
        return []
    if isinstance(val, (list, tuple)):
        items = list(val)
    else:
        s = str(val)
        for sep in ("|", ";", "\n"):
            s = s.replace(sep, ",")
        items = s.split(",")
    return [t.strip() for t in items if t and t.strip() and len(t.strip()) > 2]


def run():
    log_lines = []

    def log(msg=""):
        print(msg)
        log_lines.append(str(msg))

    log("=" * 70)
    log("DISEASE ENTITY-RESOLUTION AUDIT  -  MeSH mapping correctness")
    log(f"Run: {datetime.datetime.now().isoformat(timespec='seconds')}")
    log("=" * 70)

    conn = psycopg2.connect(**db_params())
    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        cur.execute("SELECT COUNT(*) AS n FROM mesh_diseases")
        total = cur.fetchone()["n"]
        cur.execute("SELECT mesh_id, heading, entry_terms FROM mesh_diseases "
                    "WHERE heading IS NOT NULL ORDER BY random() LIMIT %s", (SAMPLE_N,))
        sample = cur.fetchall()
    conn.close()
    log(f"\nmesh_diseases total: {total:,}   audited sample: {len(sample)}")

    def resolved_id(query):
        try:
            r = resolve_disease(query)
            return (r[0].get("mesh_id") if r else None), len(r or [])
        except Exception:
            return None, 0

    # 1. Heading round-trip
    h_ok = h_amb = 0
    for d in sample:
        rid, n = resolved_id(d["heading"])
        if rid == d["mesh_id"]:
            h_ok += 1
        if n > 1:
            h_amb += 1
    heading_pct = round(100.0 * h_ok / len(sample), 1) if sample else 0

    # 2. Synonym (entry-term) round-trip - the real test
    syn_total = syn_ok = 0
    syn_misses = []
    for d in sample:
        terms = _split_entry_terms(d.get("entry_terms"))
        if not terms:
            continue
        term = random.choice(terms)
        rid, _ = resolved_id(term)
        syn_total += 1
        if rid == d["mesh_id"]:
            syn_ok += 1
        elif len(syn_misses) < 15:
            syn_misses.append({"term": term, "expected": d["mesh_id"],
                               "expected_heading": d["heading"], "got": rid})
    syn_pct = round(100.0 * syn_ok / syn_total, 1) if syn_total else 0

    # 3. Alias resolution
    golden_ok = 0
    golden_detail = []
    for term, expect in GOLDEN_SET.items():
        try:
            r = resolve_disease(term)
            head = (r[0].get("heading") if r else "") or ""
            ok = expect.lower() in head.lower()
        except Exception:
            head, ok = "", False
        golden_ok += 1 if ok else 0
        golden_detail.append({"input": term, "expected": expect, "resolved": head, "ok": ok})
    golden_pct = round(100.0 * golden_ok / len(GOLDEN_SET), 1)

    ambiguity_pct = round(100.0 * h_amb / len(sample), 1) if sample else 0

    log(f"\n[1] Heading round-trip (consistency) : {h_ok}/{len(sample)}  {heading_pct}%")
    log(f"[2] Synonym round-trip (consistency) : {syn_ok}/{syn_total}  {syn_pct}%   (silent-killer check)")
    log(f"[3] GOLDEN SET (clinical correctness): {golden_ok}/{len(GOLDEN_SET)}  {golden_pct}%   (vs known ground truth)")
    log(f"    Ambiguous resolutions: {ambiguity_pct}% of sample returned >1 descriptor")
    if syn_misses:
        log("\n    Synonym mismatches (sample):")
        for m in syn_misses[:8]:
            log(f"      '{m['term']}' -> got {m['got']}, expected {m['expected']} ({m['expected_heading']})")

    findings = []
    if syn_pct < 90:
        findings.append({
            "id": "RES-01", "severity": "MEDIUM",
            "title": "Synonym resolution below 90%",
            "detail": f"Only {syn_pct}% of MeSH entry-term synonyms resolve to the correct descriptor. "
                      f"Mis-mapped synonyms silently route to the wrong disease.",
            "capa": "Review synonym-match precedence in resolve_disease; prefer exact entry-term hits.",
        })
    findings.append({
        "id": "RES-00", "severity": "INFO",
        "title": "Disease resolution: consistency + clinical correctness",
        "detail": f"Two distinct measures over {total:,} MeSH diseases. Consistency (round-trip identity): "
                  f"heading {heading_pct}%, synonym {syn_pct}% on a {len(sample)}-descriptor sample - proves "
                  f"the resolver is stable/deterministic. Correctness: {golden_pct}% against a hand-curated "
                  f"golden set of {len(GOLDEN_SET)} tricky synonyms/abbreviations with externally-known MeSH "
                  f"ground truth - proves names map to the CLINICALLY correct concept, not just consistently.",
        "capa": None,
    })

    result = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "audit": "MeSH disease entity-resolution (consistency + golden-set correctness)",
        "mesh_total": total,
        "sample_size": len(sample),
        "heading_roundtrip_pct": heading_pct,
        "synonym_roundtrip_pct": syn_pct,
        "synonym_pairs_tested": syn_total,
        "golden_set_accuracy_pct": golden_pct,
        "golden_set_size": len(GOLDEN_SET),
        "ambiguity_pct": ambiguity_pct,
        "golden_detail": golden_detail,
        "synonym_misses": syn_misses,
        "findings": findings,
    }
    (HERE / "resolution_results.json").write_text(json.dumps(result, indent=2))
    (HERE / "resolution_run.log").write_text("\n".join(log_lines))
    log(f"\n  Wrote: {HERE / 'resolution_results.json'}")
    return result


if __name__ == "__main__":
    try:
        run()
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(1)
