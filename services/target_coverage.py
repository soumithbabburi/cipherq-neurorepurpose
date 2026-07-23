"""
Target-coverage / mandatory-intersection gate  —  polygenic-disease scoring.
═══════════════════════════════════════════════════════════════════════════════
Step 4 of the platform architecture brief. For a complex / polygenic disease, a
drug that hits only ONE of several driver targets is rarely a viable candidate —
yet a naive target score can rank it beside a drug that covers the whole driver
set. This gate measures, generally and per-pair, how much of the disease's DRIVER
target set a drug actually covers, and applies a bounded penalty when coverage is
low AND the disease is genuinely polygenic.

Distinct from the existing "target overlap" sub-score, which is DRUG-centric
(what fraction of the DRUG's targets are disease genes). This is DISEASE-centric
(of the disease's DRIVER targets, how many does the drug hit) — the "mandatory
intersection" the brief asks for.

Key ideas, all data-driven (Open Targets genetic-association weights), no
hardcoded disease triads:
  • DRIVERS   = the disease's meaningfully-associated genes (weight ≥ 30% of the
                top gene's weight), capped.
  • POLYGENIC = effective number of drivers (inverse-Simpson of the weights). One
                dominant driver → monogenic → single-target hit suffices, NO
                penalty. Weight spread across many comparable drivers → polygenic.
  • COVERAGE  = weighted fraction of the driver set the drug hits (hitting a
                high-weight driver counts for more than a minor one).
  • PENALTY   = scales with polygenicity × (1 − coverage). Only applies to a
                MECHANISTIC candidate (hits ≥ 1 driver); a candidate surfaced by
                non-mechanistic (clinical/literature) evidence is left alone.

FAIL-SOFT: no genetic weights → multiplier 1.0 (never penalise on missing data).
"""
from __future__ import annotations

from typing import Dict, List

# A gene counts as a driver if its weight is ≥ this fraction of the top gene's.
_DRIVER_MIN_FRAC = 0.30
_DRIVER_CAP = 12
# Effective-driver-count (inverse-Simpson) ramp: <1.5 → monogenic (no penalty),
# ≥4.5 → fully polygenic.
_NEFF_MIN, _NEFF_FULL = 1.5, 4.5
_MAX_PENALTY = 0.55          # strongest coverage knock-down (at pg=1, coverage=0)
_MULTIPLIER_FLOOR = 0.30


def _effective_drivers(weights: List[float]) -> float:
    """Inverse-Simpson effective number of drivers (1 = one dominant driver)."""
    s = sum(weights)
    if s <= 0:
        return 0.0
    sq = sum(w * w for w in weights)
    return (s * s) / sq if sq > 0 else 0.0


def assess_coverage(drug_genes: List[str],
                    disease_weights: Dict[str, float]) -> Dict:
    """Assess how completely a drug covers a disease's driver target set.

    disease_weights: {GENE_UPPER: genetic_association_score 0..1}.
    Returns {multiplier, penalized, coverage, n_drivers, n_hit, effective_drivers,
             driver_genes, hit_genes, flags}. Fail-soft → multiplier 1.0.
    """
    result = {"multiplier": 1.0, "penalized": False, "coverage": None,
              "n_drivers": 0, "n_hit": 0, "effective_drivers": 0.0,
              "driver_genes": [], "hit_genes": [], "flags": []}
    if not drug_genes or not disease_weights:
        return result

    items = sorted(((g.upper(), float(w)) for g, w in disease_weights.items()
                    if w and w > 0), key=lambda kv: -kv[1])
    if not items:
        return result
    top_w = items[0][1]
    drivers = [(g, w) for g, w in items if w >= _DRIVER_MIN_FRAC * top_w][:_DRIVER_CAP]
    if len(drivers) <= 1:                      # monogenic → single target suffices
        result["n_drivers"] = len(drivers)
        result["driver_genes"] = [g for g, _ in drivers]
        return result

    drug_set = {g.upper() for g in drug_genes}
    hit = [(g, w) for g, w in drivers if g in drug_set]
    total_w = sum(w for _, w in drivers)
    hit_w = sum(w for _, w in hit)
    coverage = round(hit_w / total_w, 4) if total_w else 0.0

    result.update(
        n_drivers=len(drivers), n_hit=len(hit),
        effective_drivers=round(_effective_drivers([w for _, w in drivers]), 2),
        driver_genes=[g for g, _ in drivers], hit_genes=[g for g, _ in hit],
        coverage=coverage,
    )

    # Only a mechanistic candidate (hits ≥1 driver) is judged on coverage; a
    # candidate with no driver overlap is surfaced by other evidence and left alone.
    if not hit:
        return result

    neff = result["effective_drivers"]
    pg = max(0.0, min(1.0, (neff - _NEFF_MIN) / (_NEFF_FULL - _NEFF_MIN)))
    if pg <= 0:                                # not polygenic enough → no penalty
        return result

    multiplier = 1.0 - pg * (1.0 - coverage) * _MAX_PENALTY
    multiplier = max(_MULTIPLIER_FLOOR, round(multiplier, 3))
    result["multiplier"] = multiplier
    result["penalized"] = multiplier < 0.999
    if result["penalized"]:
        result["flags"].append(
            f"Hits {len(hit)} of {len(drivers)} disease driver targets "
            f"({round(coverage*100)}% weighted coverage); polygenic indication "
            f"(~{neff:.1f} effective drivers) — incomplete target coverage penalised.")
    return result
