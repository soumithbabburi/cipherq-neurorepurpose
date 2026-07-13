"""
Lead-viability funnel  —  potency + therapeutic window on top of plausibility.
═══════════════════════════════════════════════════════════════════════════════
The DWPC graph model finds MECHANISTICALLY PLAUSIBLE pairs but cannot tell success
from failure (validated: AUC 0.42 vs real trial failures) — efficacy/safety/PK are
not in the graph. This module adds the physical gates that ELIMINATE what cannot
work, turning a plausible hypothesis into a lab-ready candidate:

  Step 1  PLAUSIBILITY   P(treats) >= 0.15            (services/repurposing_predictor.py)
  Step 2  POTENCY        IC50 against the disease target — drop sloppy µM binders
  Step 3  WINDOW         free Cmax / IC50 — can a SAFE dose actually reach the target?

Epistemic stance: this is an ELIMINATION funnel, not a success oracle. Removing the
physically impossible (a 50 µM binder whose safe Cmax is 0.5 µM) is reliable;
predicting the possible is not. Every gate is fail-soft — missing data → no verdict
on that step, never a fabricated penalty.

Potency from the validated chembl_33 activities (pChEMBL = -log10 molar IC50/Ki/Kd).
Window from the perfusion-limited PBPK model (free Cmax at a therapeutic dose).
"""
from __future__ import annotations

import logging
import math
from functools import lru_cache
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# pChEMBL → potency class.  pChEMBL 9 = 1 nM, 7 = 100 nM, 6 = 1 µM, 5 = 10 µM.
POTENT_PCHEMBL   = 7.0     # <= 100 nM  — high-potency lead
MODERATE_PCHEMBL = 6.0     # <= 1 µM    — workable
SLOPPY_PCHEMBL   = 5.0     # >= 10 µM   — sloppy binder, eliminate
# Therapeutic window = free Cmax / IC50 (both molar). >=1 needed to engage; >=10 real margin.
WINDOW_STRONG = 10.0
WINDOW_MIN    = 1.0


def _ic50_nm(pchembl: float) -> float:
    return 10.0 ** (9.0 - pchembl)          # nM


# ── Step 2: potency (chembl_33 activities) ─────────────────────────────────────

@lru_cache(maxsize=4096)
def _best_pchembl(drug_ref: str, genes_key: str) -> Optional[tuple]:
    """(best_pchembl, gene, source) for a drug, from chembl_33, via a tiered fallback
    that trades specificity for coverage — each tier is labelled honestly:

      1. 'disease-target'   — best pChEMBL against ONE of the disease's driver genes
                              (the hypothesis target; most specific).
      2. 'mechanism-target' — best pChEMBL against the drug's ANNOTATED mechanism-of-
                              action target (its known pharmacology). Big coverage lift
                              for approved drugs whose activity isn't tagged to the OT gene.
      3. 'best-measured'    — the molecule's single most potent measured activity at ANY
                              single-protein target — the pure sloppy-binder check (a drug
                              whose BEST IC50 is >10 µM is a weak binder, period).

    drug_ref = ChEMBL id (preferred) or name. None if the drug has no measured activity."""
    genes = [g for g in genes_key.split("|") if g]
    if not drug_ref:
        return None
    try:
        import psycopg2
        from config import db_params
        p = db_params(); p["dbname"] = "chembl_33"
        conn = psycopg2.connect(**p)
    except Exception as e:
        logger.debug(f"potency: chembl_33 unavailable ({e})")
        return None
    by_id = drug_ref.upper().startswith("CHEMBL")
    where = "md.chembl_id = %s" if by_id else "md.pref_name ILIKE %s"
    try:
        with conn.cursor() as cur:
            # Tier 1 — against the disease driver genes
            if genes:
                cur.execute(
                    f"""
                    SELECT UPPER(cs.component_synonym) AS gene, MAX(act.pchembl_value) AS best
                    FROM molecule_dictionary md
                    JOIN activities act        ON act.molregno = md.molregno
                    JOIN assays a              ON a.assay_id = act.assay_id
                    JOIN target_dictionary td  ON td.tid = a.tid
                    JOIN target_components tc  ON tc.tid = td.tid
                    JOIN component_synonyms cs ON cs.component_id = tc.component_id
                    WHERE {where} AND cs.syn_type = 'GENE_SYMBOL'
                      AND UPPER(cs.component_synonym) = ANY(%s)
                      AND act.pchembl_value IS NOT NULL AND td.target_type = 'SINGLE PROTEIN'
                    GROUP BY UPPER(cs.component_synonym) ORDER BY best DESC NULLS LAST LIMIT 1
                    """,
                    (drug_ref, [g.upper() for g in genes]))
                row = cur.fetchone()
                if row and row[1] is not None:
                    return (float(row[1]), row[0], "disease-target")

            # Tier 2 — against the drug's annotated mechanism-of-action target
            cur.execute(
                f"""
                SELECT UPPER(cs.component_synonym) AS gene, MAX(act.pchembl_value) AS best
                FROM molecule_dictionary md
                JOIN drug_mechanism dm     ON dm.molregno = md.molregno
                JOIN activities act        ON act.molregno = md.molregno
                JOIN assays a              ON a.assay_id = act.assay_id AND a.tid = dm.tid
                JOIN target_components tc  ON tc.tid = dm.tid
                JOIN component_synonyms cs ON cs.component_id = tc.component_id
                WHERE {where} AND cs.syn_type = 'GENE_SYMBOL' AND act.pchembl_value IS NOT NULL
                GROUP BY UPPER(cs.component_synonym) ORDER BY best DESC NULLS LAST LIMIT 1
                """,
                (drug_ref,))
            row = cur.fetchone()
            if row and row[1] is not None:
                return (float(row[1]), row[0], "mechanism-target")

            # Tier 3 — the molecule's best CORROBORATED potency (sloppy-binder check).
            # Robust to single-assay artifacts: per single-protein target take the MEDIAN
            # pChEMBL over >=2 measurements, then the max of those target-medians. This
            # stops one spurious nM datapoint from making a genuine mM binder (metformin)
            # look potent.
            cur.execute(
                f"""
                SELECT MAX(tm) FROM (
                    SELECT PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY act.pchembl_value) AS tm
                    FROM molecule_dictionary md
                    JOIN activities act       ON act.molregno = md.molregno
                    JOIN assays a             ON a.assay_id = act.assay_id
                    JOIN target_dictionary td ON td.tid = a.tid AND td.target_type = 'SINGLE PROTEIN'
                    WHERE {where} AND act.pchembl_value IS NOT NULL
                    GROUP BY a.tid HAVING COUNT(*) >= 2
                ) t
                """,
                (drug_ref,))
            row = cur.fetchone()
            if row and row[0] is not None:
                return (float(row[0]), "", "best-measured")
    except Exception as e:
        logger.debug(f"potency query failed: {e}")
    finally:
        conn.close()
    return None


def potency(drug_ref: str, genes: List[str]) -> Dict:
    """Best measured potency of a drug against the disease-relevant target(s)."""
    out = {"covered": False, "pchembl": None, "ic50_nm": None, "target_gene": None,
           "klass": "unknown", "pass": None, "source": None,
           "note": "No measured potency (IC50/Ki) for this drug against the target(s)."}
    genes = sorted({g.strip().upper() for g in (genes or []) if g and g.strip()})
    hit = _best_pchembl(drug_ref or "", "|".join(genes))
    if not hit:
        return out
    pch, gene, source = hit
    ic50 = _ic50_nm(pch)
    if pch >= POTENT_PCHEMBL:
        klass, passed = "potent", True
    elif pch >= MODERATE_PCHEMBL:
        klass, passed = "moderate", True
    elif pch >= SLOPPY_PCHEMBL:
        klass, passed = "weak", True          # caution, not eliminated
    else:
        klass, passed = "sloppy", False       # >=10 µM — eliminate
    _srclabel = {"disease-target": "disease target", "mechanism-target": "mechanism target",
                 "best-measured": "best measured target"}.get(source, source)
    out.update(covered=True, pchembl=round(pch, 2), ic50_nm=round(ic50, 2),
               target_gene=gene or None, klass=klass, source=source, **{"pass": passed},
               note=(f"IC50 ~{_fmt_conc(ic50)}"
                     + (f" at {gene}" if gene else "") + f" (pChEMBL {pch:.1f}; {_srclabel}) — {klass} binder."
                     + ("" if passed else " Sloppy µM binder — off-target risk; eliminated.")))
    return out


def _fmt_conc(nm: float) -> str:
    if nm >= 1000:
        return f"{nm/1000:.1f} µM"
    if nm >= 1:
        return f"{nm:.0f} nM"
    return f"{nm*1000:.0f} pM"


# ── Step 3: therapeutic window (PBPK free Cmax vs IC50) ─────────────────────────

def therapeutic_window(drug_name: str, ic50_nm: float, *, mw: float = 350.0,
                       logp: float = 2.5, smiles: str = "",
                       dose_mg: float = 100.0) -> Dict:
    """Can a therapeutic dose reach the target's IC50 as FREE drug? margin = free Cmax / IC50."""
    out = {"covered": False, "free_cmax_nm": None, "margin": None, "klass": "unknown",
           "pass": None, "note": "PBPK exposure unavailable — window not assessed."}
    if not ic50_nm or ic50_nm <= 0:
        return out
    try:
        if smiles:
            try:
                from rdkit import Chem
                from rdkit.Chem import Descriptors
                m = Chem.MolFromSmiles(smiles)
                if m:
                    mw = Descriptors.MolWt(m); logp = Descriptors.MolLogP(m)
            except Exception:
                pass
        from services.pbpk_simulation import PBPKSimulator
        sim = PBPKSimulator()
        res = sim.simulate_drug_exposure(drug_name, molecular_weight=mw, logp=logp,
                                         dose_mg=dose_mg, binding_affinity=ic50_nm)
        cmax_ng_ml = float(res.get("cmax_plasma") or 0.0)
        fu = float((res.get("adme_ui") or {}).get("fu") or 1.0)
    except Exception as e:
        logger.debug(f"window PBPK failed for {drug_name}: {e}")
        return out
    if cmax_ng_ml <= 0 or mw <= 0:
        return out
    # ng/ml → nM : (ng/ml)*1000/MW ; free = ×fu
    free_cmax_nm = cmax_ng_ml * 1000.0 / mw * fu
    margin = free_cmax_nm / ic50_nm
    if margin >= WINDOW_STRONG:
        klass, passed = "strong", True
    elif margin >= WINDOW_MIN:
        klass, passed = "marginal", True
    else:
        klass, passed = "unreachable", False
    out.update(covered=True, free_cmax_nm=round(free_cmax_nm, 2), margin=round(margin, 2),
               klass=klass, **{"pass": passed},
               note=(f"Free Cmax ~{_fmt_conc(free_cmax_nm)} vs IC50 {_fmt_conc(ic50_nm)} "
                     f"→ {margin:.1f}× margin ({klass})."
                     + ("" if passed else " Safe dose cannot reach the target — eliminated.")))
    return out


# ── The funnel ─────────────────────────────────────────────────────────────────

def assess(drug_name: str, chembl_id: str, disease: str, genes: List[str], *,
           smiles: str = "", plausibility_p: Optional[float] = None,
           run_window: bool = True) -> Dict:
    """Full lead-viability funnel for a (drug, disease) hypothesis.

    Returns {steps:{plausibility,potency,window}, verdict, eliminated, reason}.
    Verdict ∈ lab-ready | needs-optimization | eliminated | insufficient-data.
    """
    drug_ref = chembl_id or drug_name
    pot = potency(drug_ref, genes)
    win = {"covered": False, "klass": "unknown", "pass": None}
    if run_window and pot.get("covered") and pot.get("ic50_nm"):
        win = therapeutic_window(drug_name, pot["ic50_nm"], smiles=smiles)

    plaus_pass = None if plausibility_p is None else (plausibility_p >= 0.15)
    # Only DISEASE-SPECIFIC or MECHANISM potency is success-relevant evidence (validated:
    # AUC 0.55 vs 0.46 for the generic best-measured tier). Generic potency can still
    # ELIMINATE a sloppy binder, but must not upgrade a pair to "lab-ready".
    _specific = pot.get("source") in ("disease-target", "mechanism-target")

    # Hard eliminations first (physics says no — valid at any potency source).
    if pot.get("pass") is False:
        verdict, elim, reason = "eliminated", True, "Sloppy µM binder (potency gate)"
    elif win.get("pass") is False:
        verdict, elim, reason = "eliminated", True, "Safe dose cannot reach the target (window gate)"
    elif not pot.get("covered"):
        verdict, elim, reason = "insufficient-data", False, "No measured potency for this pair"
    elif (pot.get("klass") == "potent" and win.get("klass") == "strong"
          and _specific and plaus_pass is not False):
        verdict, elim, reason = "lab-ready", False, "Potent at a disease/mechanism target, safe window, plausible"
    elif not _specific:
        verdict, elim, reason = ("needs-optimization", False,
                                 "Potency is generic (best-measured), not disease-specific — "
                                 "eliminates sloppy binders only")
    else:
        verdict, elim, reason = "needs-optimization", False, "Passes gates but not a clean sweep"

    return {
        "steps": {
            "plausibility": {"probability": plausibility_p,
                             "pass": plaus_pass, "cutoff": 0.15},
            "potency": pot,
            "window": win,
        },
        "verdict": verdict, "eliminated": elim, "reason": reason,
    }
