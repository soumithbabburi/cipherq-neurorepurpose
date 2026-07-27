"""
Shared Bioactivity Data Validation  -  chembl_33
════════════════════════════════════════════════
Validates the ChEMBL bioactivity/affinity data (Ki / IC50 / EC50 / Kd / ...)
that BOTH platforms read from the shared `chembl_33` database on localhost:5433:
    - CompoundIQ / POZ (screening / Cmax)        -> chembl_33
    - RepurposeIQ (docking, PBPK, repurposing)   -> chembl_33 affinities

POZ's existing `validation/validate_data.py` covers the Cmax universe only.
This script covers the AFFINITY data - the part a pharma reviewer means when
they say "ChEMBL has noise and discrepancies." It does not claim the data is
clean; it MEASURES the noise and bounds it, using ChEMBL's own quality fields.

Structured the same way as POZ's SOP-001 (ALCOA+ / GAMP 5): this script IS the
executable test protocol. Read-only - SELECT queries only, never modifies data.

ALCOA+ = Attributable, Legible, Contemporaneous, Original, Accurate,
         + Complete, Consistent, Enduring, Available

Usage:
    python validation/validate_bioactivity.py

Outputs (written next to this file):
    validation_results_bioactivity.json   - machine-readable results + scorecard
    validation_run_bioactivity.log        - console transcript
"""

import os
import sys
import json
import datetime
from pathlib import Path

from sqlalchemy import create_engine, text
from sqlalchemy.engine import URL

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

# Affinity endpoints we validate, each kept in its OWN lane (never mixed).
# These are -log10(molar) potency measures expressed as pchembl_value, so a
# spread of 1.0 in pchembl == 1 order of magnitude (1 log unit).
AFFINITY_TYPES = ["Ki", "IC50", "EC50", "Kd"]

# ChEMBL assay->target confidence: >=8 means a single, directly-assigned
# protein target (the gold tier). Below that the target assignment is fuzzier.
HIGH_CONFIDENCE = 8

# Published experimental reproducibility floor for public affinity data
# (~0.5 pKi/pIC50 log units; Kramer et al. 2012, J Med Chem). We compare our
# MEASURED intra-pair spread against this reference.
LITERATURE_NOISE_FLOOR_LOG = 0.5


def load_env():
    env = {}
    envfile = ROOT / ".env"
    if envfile.exists():
        for line in envfile.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            env[k.strip()] = v.strip()
    return env


def get_engine():
    """Connect to the SHARED chembl_33 DB. Creds from .env / env; never hardcoded.
    dbname is forced to chembl_33 (cipherq's .env DB_NAME points at neurorepurpose)."""
    env = load_env()

    def pick(*keys, default=""):
        for k in keys:
            v = env.get(k) or os.environ.get(k)
            if v:
                return v
        return default

    user = pick("DB_USER", "CHEMBL_DB_USER", default="babburisoumith")
    pw = pick("DB_PASSWORD", "CHEMBL_DB_PASSWORD", "PGPASSWORD")
    host = pick("DB_HOST", "CHEMBL_DB_HOST", default="localhost")
    port = pick("DB_PORT", "CHEMBL_DB_PORT", default="5433")
    # URL.create handles special characters in the password (e.g. '@') safely.
    uri = URL.create(
        "postgresql+psycopg2",
        username=user, password=pw or None,
        host=host, port=int(port), database="chembl_33",
    )
    return create_engine(uri, connect_args={"connect_timeout": 30})


def run():
    log_lines = []

    def log(msg=""):
        print(msg)
        log_lines.append(str(msg))

    eng = get_engine()

    def scalar(sql, **p):
        with eng.connect() as c:
            return c.execute(text(sql), p).scalar()

    def rows(sql, **p):
        with eng.connect() as c:
            return c.execute(text(sql), p).fetchall()

    results = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "database": "chembl_33",
        "scope": "Bioactivity/affinity data (Ki/IC50/EC50/Kd) for RepurposeIQ",
        "standard": "ALCOA+ (FDA/WHO/EU Annex 11), GAMP 5 structure",
        "reference_noise_floor_log": LITERATURE_NOISE_FLOOR_LOG,
        "checks": {},
        "findings": [],
        "scorecard": {},
    }

    def finding(fid, severity, title, detail, capa=None):
        results["findings"].append({
            "id": fid, "severity": severity, "title": title,
            "detail": detail, "capa": capa,
        })

    log("=" * 72)
    log("SHARED BIOACTIVITY DATA VALIDATION  -  chembl_33  -  ALCOA+ / GAMP 5")
    log(f"Run: {results['run_at']}")
    log("=" * 72)

    # ════════════════════════════════════════════════════════════════
    # 0. AVAILABLE / ENDURING + TIMELINESS (version stamp)
    # ════════════════════════════════════════════════════════════════
    log("\n[0] AVAILABLE / ENDURING / TIMELINESS")
    total_act = scalar("SELECT COUNT(*) FROM activities")
    ver = rows("SELECT name, creation_date FROM version LIMIT 1")
    ver_name = ver[0][0] if ver else "unknown"
    ver_date = str(ver[0][1]) if ver else "unknown"
    results["checks"]["available"] = {
        "activities_rows": total_act, "chembl_version": ver_name,
        "release_date": ver_date, "pass": total_act > 0,
    }
    log(f"    activities rows : {total_act:,}")
    log(f"    ChEMBL version  : {ver_name}  (released {ver_date})")

    # ════════════════════════════════════════════════════════════════
    # 1. SCOPE - the affinity universe, by endpoint (kept SEPARATE)
    # ════════════════════════════════════════════════════════════════
    log("\n[1] SCOPE - affinity universe by endpoint (assay types NOT mixed)")
    universe = {}
    for t in AFFINITY_TYPES:
        n_all = scalar("SELECT COUNT(*) FROM activities WHERE standard_type=:t", t=t)
        n_pc = scalar("SELECT COUNT(*) FROM activities WHERE standard_type=:t AND pchembl_value IS NOT NULL", t=t)
        universe[t] = {"rows": n_all, "with_pchembl": n_pc}
        log(f"    {t:5s}: {n_all:>9,} rows   {n_pc:>9,} with pchembl_value")
    results["checks"]["universe_by_endpoint"] = universe

    # ════════════════════════════════════════════════════════════════
    # 2. ATTRIBUTABLE - every measurement traces to an assay + target
    # ════════════════════════════════════════════════════════════════
    log("\n[2] ATTRIBUTABLE - measurement -> assay -> target traceability")
    no_assay = scalar("""SELECT COUNT(*) FROM activities a
        WHERE a.standard_type = ANY(:types) AND a.assay_id IS NULL""", types=AFFINITY_TYPES)
    no_target = scalar("""SELECT COUNT(*) FROM activities a
        JOIN assays ass ON a.assay_id = ass.assay_id
        WHERE a.standard_type = ANY(:types) AND ass.tid IS NULL""", types=AFFINITY_TYPES)
    results["checks"]["attributable"] = {
        "activities_without_assay": no_assay,
        "activities_without_target": no_target,
        "pass": (no_assay == 0),
    }
    log(f"    affinity rows with no assay link : {no_assay:,}")
    log(f"    affinity rows with no target     : {no_target:,}")
    if no_target > 0:
        finding("ATTR-B01", "INFO", "Affinity measurements with no assigned target",
                f"{no_target:,} affinity rows link to an assay that has no tid (target). "
                f"These cannot anchor a compound-target pair and are excluded from pair analysis.",
                "Exclude from target-based scoring; documented, not silently dropped.")

    # ════════════════════════════════════════════════════════════════
    # 3. VALIDITY - use ChEMBL's OWN quality flags
    # ════════════════════════════════════════════════════════════════
    log("\n[3] VALIDITY - ChEMBL native quality flags (confidence / validity / censoring)")
    validity = {}
    for t in AFFINITY_TYPES:
        base = "FROM activities a JOIN assays ass ON a.assay_id=ass.assay_id WHERE a.standard_type=:t AND a.pchembl_value IS NOT NULL"
        tot = scalar(f"SELECT COUNT(*) {base}", t=t)
        hi = scalar(f"SELECT COUNT(*) {base} AND ass.confidence_score >= :c", t=t, c=HIGH_CONFIDENCE)
        flagged = scalar(f"SELECT COUNT(*) {base} AND a.data_validity_comment IS NOT NULL", t=t)
        censored = scalar(f"SELECT COUNT(*) {base} AND a.standard_relation IS DISTINCT FROM '='", t=t)
        dupe = scalar(f"SELECT COUNT(*) {base} AND a.potential_duplicate = 1", t=t)
        oor = scalar(f"SELECT COUNT(*) {base} AND (a.pchembl_value <= 0 OR a.pchembl_value > 14)", t=t)
        validity[t] = {
            "total_with_pchembl": tot,
            "high_confidence_ge8": hi,
            "high_confidence_pct": round(100.0 * hi / tot, 1) if tot else 0,
            "data_validity_flagged": flagged,
            "censored_relation": censored,
            "potential_duplicate": dupe,
            "pchembl_out_of_range": oor,
        }
        log(f"    {t:5s}: hi-conf {hi:>8,} ({validity[t]['high_confidence_pct']:>5}%)  "
            f"validity-flagged {flagged:>6,}  censored {censored:>7,}  dup {dupe:>6,}  oor {oor}")
    results["checks"]["validity"] = validity
    tot_flagged = sum(v["data_validity_flagged"] for v in validity.values())
    if tot_flagged > 0:
        finding("VAL-B01", "MEDIUM", "ChEMBL-flagged suspect measurements present",
                f"{tot_flagged:,} affinity rows carry a non-null data_validity_comment "
                f"(ChEMBL's own 'this value looks wrong' flag). Excluded from the certified tier.",
                "Filter on data_validity_comment IS NULL for certified affinities; keep raw rows flagged.")

    # ════════════════════════════════════════════════════════════════
    # 4. NOISE / REPRODUCIBILITY - the core 'how noisy is ChEMBL' metric
    # ════════════════════════════════════════════════════════════════
    log("\n[4] REPRODUCIBILITY - measured intra-pair spread (the noise number)")
    log(f"    (high-confidence cs>={HIGH_CONFIDENCE}, validity-clean, relation '=', >=2 measurements/pair)")
    noise = {}
    NOISE_SQL = """
    WITH clean AS (
        SELECT a.molregno, ass.tid, a.pchembl_value AS pv
        FROM activities a JOIN assays ass ON a.assay_id=ass.assay_id
        WHERE a.standard_type=:t AND a.pchembl_value IS NOT NULL
          AND ass.confidence_score >= :c
          AND a.data_validity_comment IS NULL
          AND a.standard_relation = '='
    ), g AS (
        SELECT molregno, tid, COUNT(*) n, stddev_pop(pv) sd
        FROM clean GROUP BY molregno, tid HAVING COUNT(*) >= 2
    )
    SELECT COUNT(*) pairs,
           round(avg(sd)::numeric,3) mean_sd,
           round((percentile_cont(0.5) WITHIN GROUP (ORDER BY sd))::numeric,3) median_sd,
           round((percentile_cont(0.95) WITHIN GROUP (ORDER BY sd))::numeric,3) p95_sd,
           round((100.0*sum(CASE WHEN sd>1 THEN 1 ELSE 0 END)/NULLIF(COUNT(*),0))::numeric,2) pct_gt1log
    FROM g
    """
    for t in AFFINITY_TYPES:
        r = rows(NOISE_SQL, t=t, c=HIGH_CONFIDENCE)
        if r and r[0][0]:
            pairs, mean_sd, median_sd, p95_sd, pct = r[0]
            noise[t] = {
                "multi_measured_pairs": int(pairs),
                "mean_log_sd": float(mean_sd) if mean_sd is not None else None,
                "median_log_sd": float(median_sd) if median_sd is not None else None,
                "p95_log_sd": float(p95_sd) if p95_sd is not None else None,
                "pct_pairs_disagree_gt_1log": float(pct) if pct is not None else None,
                "vs_literature_floor": "within" if (mean_sd is not None and float(mean_sd) <= LITERATURE_NOISE_FLOOR_LOG) else "above",
            }
            log(f"    {t:5s}: pairs {int(pairs):>7,}  median {float(median_sd):.3f}  "
                f"mean {float(mean_sd):.3f}  p95 {float(p95_sd):.3f} log  "
                f">1log {float(pct):.2f}%  ({noise[t]['vs_literature_floor']} {LITERATURE_NOISE_FLOOR_LOG} floor)")
        else:
            noise[t] = {"multi_measured_pairs": 0}
            log(f"    {t:5s}: no multi-measured high-confidence pairs")
    results["checks"]["reproducibility"] = noise

    # sample the worst-disagreeing Ki pairs (named) for the flagged list
    log("\n    Worst-disagreeing high-confidence Ki pairs (sample, flagged for review):")
    worst = rows("""
        WITH clean AS (
            SELECT a.molregno, ass.tid, a.pchembl_value pv
            FROM activities a JOIN assays ass ON a.assay_id=ass.assay_id
            WHERE a.standard_type='Ki' AND a.pchembl_value IS NOT NULL
              AND ass.confidence_score >= :c AND a.data_validity_comment IS NULL
              AND a.standard_relation='='
        ), g AS (
            SELECT molregno, tid, COUNT(*) n, stddev_pop(pv) sd,
                   min(pv) lo, max(pv) hi
            FROM clean GROUP BY molregno, tid HAVING COUNT(*) >= 3
        )
        SELECT md.chembl_id, COALESCE(md.pref_name,'(unnamed)') name,
               td.pref_name target, g.n, round(g.sd::numeric,2) sd,
               round(g.lo::numeric,2) lo, round(g.hi::numeric,2) hi
        FROM g JOIN molecule_dictionary md ON g.molregno=md.molregno
               JOIN target_dictionary td ON g.tid=td.tid
        ORDER BY g.sd DESC LIMIT 15
    """, c=HIGH_CONFIDENCE)
    flagged_pairs = []
    for r in worst:
        flagged_pairs.append({
            "chembl_id": r[0], "drug": r[1], "target": r[2],
            "n_measurements": int(r[3]), "log_sd": float(r[4]),
            "pchembl_min": float(r[5]), "pchembl_max": float(r[6]),
        })
        log(f"        {r[0]:14s} {str(r[1])[:22]:22s} {str(r[2])[:26]:26s} "
            f"n={r[3]} sd={float(r[4]):.2f} range {float(r[5]):.2f}-{float(r[6]):.2f}")
    results["checks"]["flagged_disagreeing_pairs_ki"] = flagged_pairs
    if flagged_pairs:
        finding("NOISE-B01", "INFO", "High-disagreement compound-target pairs identified",
                "Compound-target pairs whose replicate affinity measurements disagree by "
                ">1 log are flagged. They are NOT deleted; the certified value is the median "
                "and the spread is carried as the uncertainty band on that pair.",
                "Surface median +/- spread (not a point value) for these pairs; show N and range.")

    # ════════════════════════════════════════════════════════════════
    # 5. UNIQUENESS - duplicate accounting + aggregation plan
    # ════════════════════════════════════════════════════════════════
    log("\n[5] UNIQUENESS - duplicate measurements -> aggregate to median (keep raw)")
    uniq = {}
    for t in AFFINITY_TYPES:
        r = rows("""
            WITH clean AS (
                SELECT a.molregno, ass.tid
                FROM activities a JOIN assays ass ON a.assay_id=ass.assay_id
                WHERE a.standard_type=:t AND a.pchembl_value IS NOT NULL
                  AND ass.confidence_score >= :c AND a.data_validity_comment IS NULL
                  AND a.standard_relation='='
            )
            SELECT COUNT(*) measurements, COUNT(DISTINCT (molregno,tid)) pairs
            FROM clean
        """, t=t, c=HIGH_CONFIDENCE)
        m, p = (int(r[0][0]), int(r[0][1])) if r and r[0][0] else (0, 0)
        uniq[t] = {"measurements": m, "distinct_pairs": p,
                   "avg_measurements_per_pair": round(m / p, 2) if p else 0}
        log(f"    {t:5s}: {m:>8,} measurements over {p:>8,} pairs "
            f"({uniq[t]['avg_measurements_per_pair']}x avg)")
    results["checks"]["uniqueness"] = uniq

    # ════════════════════════════════════════════════════════════════
    # 6. COMPLETENESS - field coverage on the affinity data
    # ════════════════════════════════════════════════════════════════
    log("\n[6] COMPLETENESS - required-field coverage (all affinity rows)")
    base_all = "FROM activities a JOIN assays ass ON a.assay_id=ass.assay_id WHERE a.standard_type = ANY(:types)"
    tot_aff = scalar(f"SELECT COUNT(*) {base_all}", types=AFFINITY_TYPES)
    has_pc = scalar(f"SELECT COUNT(*) {base_all} AND a.pchembl_value IS NOT NULL", types=AFFINITY_TYPES)
    has_units = scalar(f"SELECT COUNT(*) {base_all} AND a.standard_units IS NOT NULL", types=AFFINITY_TYPES)
    has_cs = scalar(f"SELECT COUNT(*) {base_all} AND ass.confidence_score IS NOT NULL", types=AFFINITY_TYPES)
    comp = {
        "total_affinity_rows": tot_aff,
        "pchembl_coverage_pct": round(100.0 * has_pc / tot_aff, 1) if tot_aff else 0,
        "units_coverage_pct": round(100.0 * has_units / tot_aff, 1) if tot_aff else 0,
        "confidence_coverage_pct": round(100.0 * has_cs / tot_aff, 1) if tot_aff else 0,
    }
    results["checks"]["completeness"] = comp
    log(f"    total affinity rows   : {tot_aff:,}")
    log(f"    pchembl_value coverage: {comp['pchembl_coverage_pct']}%")
    log(f"    units coverage        : {comp['units_coverage_pct']}%")
    log(f"    confidence coverage   : {comp['confidence_coverage_pct']}%")

    # ════════════════════════════════════════════════════════════════
    # SCORECARD - six dimensions for the chembl_33 bioactivity dataset
    # ════════════════════════════════════════════════════════════════
    log("\n" + "=" * 72)
    log("DATA-QUALITY SCORECARD  -  chembl_33 bioactivity")
    log("=" * 72)
    hc_overall = (sum(v["high_confidence_ge8"] for v in validity.values())
                  / max(1, sum(v["total_with_pchembl"] for v in validity.values())))
    worst_mean_noise = max((n.get("mean_log_sd") or 0) for n in noise.values()) if noise else None
    scorecard = {
        "completeness": {"metric": "pchembl coverage", "value_pct": comp["pchembl_coverage_pct"]},
        "validity": {"metric": "high-confidence (cs>=8) share", "value_pct": round(100 * hc_overall, 1)},
        "accuracy": {"metric": "measured noise vs 0.5 log floor",
                     "worst_mean_log_sd": worst_mean_noise,
                     "status": "within floor" if (worst_mean_noise is not None and worst_mean_noise <= LITERATURE_NOISE_FLOOR_LOG) else "review"},
        "consistency": {"metric": "assay types kept separate", "status": "enforced (per-endpoint lanes)"},
        "uniqueness": {"metric": "duplicates aggregated to median", "status": "planned (raw retained)"},
        "timeliness": {"metric": "source version", "value": f"{ver_name} ({ver_date})"},
    }
    results["scorecard"] = scorecard
    for dim, d in scorecard.items():
        log(f"    {dim:13s}: {json.dumps(d)}")

    # ──────────────────────────────────────────────────────────────
    log("\n" + "=" * 72)
    log("SUMMARY OF FINDINGS")
    log("=" * 72)
    for f in results["findings"]:
        log(f"  [{f['severity']:6s}] {f['id']}: {f['title']}")
    log(f"\n  Total findings: {len(results['findings'])}")

    (HERE / "validation_results_bioactivity.json").write_text(json.dumps(results, indent=2))
    (HERE / "validation_run_bioactivity.log").write_text("\n".join(log_lines))
    log(f"\n  Wrote: {HERE / 'validation_results_bioactivity.json'}")
    log(f"  Wrote: {HERE / 'validation_run_bioactivity.log'}")
    return results


if __name__ == "__main__":
    try:
        run()
    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit(1)
