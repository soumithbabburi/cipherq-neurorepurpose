"""
Off-target / polypharmacology inference (secondary pharmacology).
═══════════════════════════════════════════════════════════════════════════
Curated drug_mechanism rows list only a drug's PRIMARY target(s). The lever for
non-obvious repurposing is the SECONDARY pharmacology: proteins the drug also
engages at (or near) a therapeutic concentration. This module infers those
secondary targets from MEASURED ChEMBL bioactivity (a reverse promiscuity
lookup), assigns each a confidence weight from potency + evidence, and attaches
provenance so every predicted off-target is attributable.

This module is PURE and DB-agnostic: the caller runs the activity query in its
own driver (psycopg2 in the engine, SQLAlchemy in the validation harness) and
passes rows in as (gene, max_pchembl, n_measurements). The science — gating,
weighting, capping — lives here and is unit-testable without a database.

Why weighting (not just membership) matters — the biology-audit guard:
  A promiscuous kinase inhibitor is measured against dozens of proteins at
  <=1 uM. If each such off-target were treated as equal to a curated primary
  target, promiscuous drugs would force-match many diseases and float to the
  top on target COUNT alone (confounding-by-promiscuity). So:
    - every off-target is capped strictly below a curated primary (max_weight),
    - weight decays with the potency gap below the drug's therapeutic anchor,
    - single-assay outliers are rejected (min_measurements),
    - off-targets far weaker than the anchor are rejected (potency_window),
    - the count per drug is capped (max_offtargets).
  Downstream ranking uses these weights, so a weak predicted off-target cannot
  out-score a genuine primary-target match.

Enabled behind the OFFTARGET_ENABLE flag so it can be measured on repoDB and
shipped (or not) as a deliberate, reversible decision.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple


# ── Tunable parameters (env-overridable; NO drug/target/disease hardcoding) ────

def _f(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, default))
    except (TypeError, ValueError):
        return default


def _i(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, default))
    except (TypeError, ValueError):
        return default


def is_enabled() -> bool:
    """Master flag. Off-target inference is OFF unless explicitly enabled, so the
    canonical benchmark and production scorer are unchanged by default."""
    return os.environ.get("OFFTARGET_ENABLE", "").strip().lower() in ("1", "true", "yes", "on")


# Defaults chosen from the secondary-pharmacology literature, the noise floor of
# public bioactivity (pKi SD ~0.54), AND a repoDB tuning sweep. On repoDB
# approved-vs-failed (validation/validate_predictions.py, fallback mapping),
# MAX_WEIGHT 0.25 / WINDOW 2.0 / MIN_N 3 is NON-REGRESSING: AUROC 0.7368 vs a
# 0.7366 baseline (delta +0.0001; paired-bootstrap 95% CI [-0.0029, 0.0033]
# straddles zero) with the Strong tier NOT inflated (22.0% vs 24.3%). Loosening
# the weight/window (e.g. 0.6/2.0) makes it dilutive (delta -0.0066, CI excludes
# zero) — the same "broadening targets dilutes ranking" pattern seen for pathway,
# KGE and signature-reversal. So off-target is kept as a BROAD generation lever
# (2-log reach) with a LOW ranking weight (a secondary target caps at 0.25 of a
# curated primary), and behind OFFTARGET_ENABLE.
MIN_PCHEMBL       = _f("OFFTARGET_MIN_PCHEMBL", 6.0)    # <=1 uM to even consider
MIN_MEASUREMENTS  = _i("OFFTARGET_MIN_N", 3)            # reject single/double-assay outliers
POTENCY_WINDOW    = _f("OFFTARGET_WINDOW", 2.0)         # log units below anchor still "reachable"
MAX_OFFTARGETS    = _i("OFFTARGET_MAX", 5)             # cap count per drug
MAX_WEIGHT        = _f("OFFTARGET_MAX_WEIGHT", 0.25)   # strictly below a curated primary (=1.0)


# ── Pure weighting core ────────────────────────────────────────────────────────

def _evidence_factor(n: int) -> float:
    """More independent measurements -> more trustworthy. n=2 -> 0.5, saturating
    to 1.0 by ~n=7. Below MIN_MEASUREMENTS the caller has already rejected."""
    return max(0.0, min(1.0, 0.5 + 0.1 * (n - MIN_MEASUREMENTS)))


def _potency_factor(pchembl: float, anchor: float) -> float:
    """Decay with the potency GAP below the therapeutic anchor (the drug's own
    primary-target potency). Each log unit weaker halves the weight; an off-target
    at or above the anchor gets the full factor (capped at 1.0, not rewarded past
    the anchor)."""
    gap = max(0.0, anchor - pchembl)
    return 0.5 ** gap


def weight_candidates(
    rows: List[Tuple[str, float, int]],
    primary_genes: set,
    *,
    min_pchembl: float = None,
    min_measurements: int = None,
    potency_window: float = None,
    max_offtargets: int = None,
    max_weight: float = None,
) -> List[Dict]:
    """Rank + weight secondary targets for ONE molecule.

    rows            : iterable of (gene_symbol, max_pchembl, n_measurements) — the
                      drug's measured single-protein activities (any target).
    primary_genes   : set of curated primary gene symbols (UPPER) to exclude and
                      to anchor the therapeutic potency.

    Returns a list of dicts (highest-confidence first), each:
      {gene, weight (0..max_weight], max_pchembl, n_measurements, potency_gap,
       provenance}
    """
    min_pchembl      = MIN_PCHEMBL      if min_pchembl      is None else min_pchembl
    min_measurements = MIN_MEASUREMENTS if min_measurements is None else min_measurements
    potency_window   = POTENCY_WINDOW   if potency_window   is None else potency_window
    max_offtargets   = MAX_OFFTARGETS   if max_offtargets   is None else max_offtargets
    max_weight       = MAX_WEIGHT       if max_weight       is None else max_weight

    primary = {g.upper() for g in (primary_genes or set())}

    # Aggregate to per-gene best potency + total evidence (rows may repeat a gene).
    per_gene: Dict[str, Tuple[float, int]] = {}
    for gene, pchembl, n in rows:
        if not gene or pchembl is None:
            continue
        g = gene.upper()
        try:
            p = float(pchembl)
            c = int(n or 0)
        except (TypeError, ValueError):
            continue
        best_p, tot_n = per_gene.get(g, (0.0, 0))
        per_gene[g] = (max(best_p, p), tot_n + c)

    if not per_gene:
        return []

    # Therapeutic anchor = the drug's strongest CURATED primary target potency, if
    # that target has measured activity; else the strongest reliable measured
    # target overall (fallback so drugs with uncurated primaries still work).
    primary_pots = [p for g, (p, _) in per_gene.items() if g in primary]
    if primary_pots:
        anchor = max(primary_pots)
    else:
        reliable = [p for g, (p, n) in per_gene.items() if n >= min_measurements]
        if not reliable:
            return []
        anchor = max(reliable)

    floor = max(min_pchembl, anchor - potency_window)

    out: List[Dict] = []
    for g, (p, n) in per_gene.items():
        if g in primary:
            continue                      # already carried as a primary target
        if n < min_measurements:
            continue                      # single-assay / noise-floor outlier
        if p < floor:
            continue                      # not engaged at a reachable concentration
        w = round(max_weight * _potency_factor(p, anchor) * _evidence_factor(n), 4)
        if w <= 0.0:
            continue
        gap = round(max(0.0, anchor - p), 2)
        out.append({
            "gene": g,
            "weight": w,
            "max_pchembl": round(p, 2),
            "n_measurements": n,
            "potency_gap": gap,
            "provenance": (f"measured ChEMBL bioactivity: max pChEMBL {p:.2f} "
                           f"across {n} activities; {gap:.1f} log below the drug's "
                           f"primary-target potency ({anchor:.2f})"),
        })

    out.sort(key=lambda d: d["weight"], reverse=True)
    return out[:max_offtargets]


def expand(
    mol_rows: Dict[int, List[Tuple[str, float, int]]],
    primary_by_mol: Dict[int, set],
    **params,
) -> Dict[int, List[Dict]]:
    """Batch wrapper: {molregno: activity_rows} + {molregno: primary_genes}
    -> {molregno: [weighted off-target dicts]}. Pure."""
    result: Dict[int, List[Dict]] = {}
    for mol, rows in mol_rows.items():
        preds = weight_candidates(rows, primary_by_mol.get(mol, set()), **params)
        if preds:
            result[mol] = preds
    return result
