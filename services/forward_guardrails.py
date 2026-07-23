"""
Forward-path development-reality guardrails (disease → drugs discovery).
════════════════════════════════════════════════════════════════════════════════
The forward discovery screen generates candidates by target/network match, which on a
saturated target rewards exactly the drugs that are ALREADY the standard of care —
e.g. asthma → ADRB2 → a list of beta 2 agonists, several of them already approved for
asthma. Target matching alone cannot tell a novel repurposing lead from a me-too of the
existing frontline drug.

These guardrails add the missing development reality, all GROUNDED in the local
chembl-derived tables (no live calls, id-independent, matched by drug name):

  • approved_for_indication — is the molecule already approved (phase >= 4) for the
    searched disease? The single true novelty disqualifier (exclude, do not rank).
  • mechanism_crowding — is the molecule's mechanism (its targets) ALREADY occupied by
    an approved drug for this disease? If every disease-relevant target it hits already
    has an approved drug for the indication, it is a me-too offering no differentiation
    (demote + flag). A candidate that reaches the disease through a target with NO
    approved drug is genuine whitespace (surfaced).

Fail-soft: a DB hiccup returns a neutral verdict, never a false disqualification.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional, Set

logger = logging.getLogger(__name__)

_MULT_CROWDED = 0.45          # every relevant target already has an approved drug: me-too
_MULT_PARTIAL = 0.75          # crowded on some targets but also hits fresh mechanism

# Non-therapeutic genes: drug-metabolising enzymes, efflux transporters, and carrier
# proteins are PK machinery hit by many drugs, not the therapeutic mechanism — they
# would inflate the "occupied" space and mislabel a novel drug as a me-too. Excluded
# from crowding on both sides.
_NON_THERAPEUTIC_PREFIX = ("CYP", "ABC", "SLCO", "UGT", "SULT")
_NON_THERAPEUTIC_EXACT = {"ALB", "AHR", "NR1I2", "NR1I3"}  # albumin, PXR/CAR/AhR (ADME regulators)


def _therapeutic(genes: Set[str]) -> Set[str]:
    return {g for g in genes
            if not g.startswith(_NON_THERAPEUTIC_PREFIX) and g not in _NON_THERAPEUTIC_EXACT}


def _q(sql: str, params=()):
    """Run a read query against the local compound DB. Returns [] on any failure."""
    try:
        from config import db_params
        import psycopg2, psycopg2.extras
        cn = psycopg2.connect(**db_params())
        cn.autocommit = True
        c = cn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
        c.execute(sql, params)
        rows = c.fetchall()
        cn.close()
        return rows
    except Exception as e:
        logger.debug("forward_guardrails query failed: %s", e)
        return []


def _disease_match(disease: str) -> str:
    """A LIKE pattern for the disease name (loose, lower-cased)."""
    return "%" + (disease or "").strip().lower() + "%"


def approved_for_indication(name: str, disease: str) -> Dict:
    """Is `name` approved (phase >= 4) for `disease` in the local indication data?
    Returns {approved, matched_disease, phase}. Grounded, name-keyed, id-independent."""
    if not name or not disease:
        return {"approved": False}
    rows = _q(
        """
        SELECT i.disease, MAX(i.max_phase) AS phase
        FROM indications i JOIN compounds co ON co.id = i.compound_id
        WHERE lower(co.name) = lower(%s) AND lower(i.disease) LIKE %s
        GROUP BY i.disease
        ORDER BY phase DESC LIMIT 1
        """,
        (name, _disease_match(disease)),
    )
    if rows and float(rows[0].get("phase") or 0) >= 4:
        return {"approved": True, "matched_disease": rows[0]["disease"],
                "phase": float(rows[0]["phase"])}
    return {"approved": False}


def drug_targets(name: str) -> Set[str]:
    """Gene symbols the drug acts on (local mechanisms table)."""
    rows = _q(
        """
        SELECT DISTINCT t.gene_symbol
        FROM mechanisms m JOIN compounds co ON co.id = m.compound_id
        JOIN targets t ON t.id = m.target_id
        WHERE lower(co.name) = lower(%s) AND t.gene_symbol IS NOT NULL
        """,
        (name,),
    )
    return {(r["gene_symbol"] or "").strip().upper() for r in rows if r.get("gene_symbol")}


def occupied_targets(disease: str) -> Set[str]:
    """Gene symbols that ALREADY have an approved (phase >= 4) drug for this disease —
    the mechanism space that is already solved for the indication."""
    rows = _q(
        """
        SELECT DISTINCT t.gene_symbol
        FROM indications i
        JOIN compounds co ON co.id = i.compound_id
        JOIN mechanisms m ON m.compound_id = co.id
        JOIN targets t ON t.id = m.target_id
        WHERE lower(i.disease) LIKE %s AND i.max_phase >= 4 AND t.gene_symbol IS NOT NULL
        """,
        (_disease_match(disease),),
    )
    return {(r["gene_symbol"] or "").strip().upper() for r in rows if r.get("gene_symbol")}


def mechanism_crowding(name: str, disease: str,
                       occupied: Optional[Set[str]] = None,
                       candidate_targets: Optional[Set[str]] = None) -> Dict:
    """Is this candidate a me-too of the approved mechanism for the disease?

    Returns {crowded, novel, multiplier, crowded_targets, novel_targets, note}.
    crowded  = every disease-relevant target it hits already has an approved drug.
    novel    = it reaches the disease through at least one un-drugged (whitespace) target.
    """
    occ = _therapeutic(occupied if occupied is not None else occupied_targets(disease))
    tgts = _therapeutic(candidate_targets if candidate_targets is not None else drug_targets(name))
    if not tgts or not occ:
        return {"crowded": False, "novel": False, "multiplier": 1.0,
                "crowded_targets": [], "novel_targets": [], "note": ""}

    crowded_t = sorted(tgts & occ)
    novel_t = sorted(tgts - occ)
    # Relevant = targets that are either occupied (shared with approved space) or the
    # drug's mechanism at large; whitespace is any target with no approved drug here.
    if crowded_t and not novel_t:
        return {"crowded": True, "novel": False, "multiplier": _MULT_CROWDED,
                "crowded_targets": crowded_t, "novel_targets": [],
                "note": ("Me-too mechanism: every target this drug hits for the indication "
                         "(" + ", ".join(crowded_t) + ") already has an approved drug for "
                         "this disease, so it offers no mechanistic differentiation.")}
    if crowded_t and novel_t:
        return {"crowded": False, "novel": True, "multiplier": _MULT_PARTIAL,
                "crowded_targets": crowded_t, "novel_targets": novel_t,
                "note": ("Partly crowded: shares " + ", ".join(crowded_t) + " with approved "
                         "drugs but also reaches " + ", ".join(novel_t) + ", which is less "
                         "exploited for this disease.")}
    return {"crowded": False, "novel": True, "multiplier": 1.0,
            "crowded_targets": [], "novel_targets": novel_t,
            "note": ("Novel mechanism for this disease: reaches it through "
                     + ", ".join(novel_t) + ", which has no approved drug for the "
                     "indication (mechanistic whitespace).")}


# ── Batched helpers (one query for a whole candidate list) ─────────────────────
def approved_names_for_indication(names: List[str], disease: str) -> Set[str]:
    """Lower-cased names in `names` that are approved (phase >= 4) for `disease` locally."""
    names = [n for n in {(_x or "").strip().lower() for _x in names} if n]
    if not names:
        return set()
    rows = _q(
        """
        SELECT DISTINCT lower(co.name) AS nm
        FROM indications i JOIN compounds co ON co.id = i.compound_id
        WHERE lower(co.name) = ANY(%s) AND lower(i.disease) LIKE %s AND i.max_phase >= 4
        """,
        (names, _disease_match(disease)),
    )
    return {r["nm"] for r in rows if r.get("nm")}


def targets_for_names(names: List[str]) -> Dict[str, Set[str]]:
    """{lower_name: {gene symbols}} for a whole candidate list in one query."""
    names = [n for n in {(_x or "").strip().lower() for _x in names} if n]
    out: Dict[str, Set[str]] = {n: set() for n in names}
    if not names:
        return out
    rows = _q(
        """
        SELECT lower(co.name) AS nm, t.gene_symbol AS g
        FROM mechanisms m JOIN compounds co ON co.id = m.compound_id
        JOIN targets t ON t.id = m.target_id
        WHERE lower(co.name) = ANY(%s) AND t.gene_symbol IS NOT NULL
        """,
        (names,),
    )
    for r in rows:
        out.setdefault(r["nm"], set()).add((r["g"] or "").strip().upper())
    return out


def apply(candidates: List[Dict], disease: str) -> Dict:
    """Apply the forward-path development-reality guardrails to a scored candidate list.

    Annotates each candidate in place with `approved_here`, `market_status`,
    `mechanism_crowding`, and a `reality_multiplier`, then partitions:
      leads   — viable, ranked by score × reality_multiplier (me-too / soft-status demoted)
      removed — approved for this indication, or withdrawn for safety (with the reason)

    Returns {"leads": [...], "removed": [...]}.
    """
    from services import market_status
    if not candidates:
        return {"leads": [], "removed": []}

    names = [c.get("name", "") for c in candidates]
    approved = approved_names_for_indication(names, disease)
    tmap = targets_for_names(names)
    occ = occupied_targets(disease)

    leads, removed = [], []
    for c in candidates:
        nm = (c.get("name", "") or "").strip()
        nml = nm.lower()

        # 1. Already approved for THIS indication → not a repurposing lead. Remove.
        if nml in approved:
            c["approved_here"] = True
            c["removed_reason"] = "Already approved for this indication"
            removed.append(c)
            continue
        c["approved_here"] = False

        # 2. Market/regulatory status. Withdrawn (safety) → remove; soft → demote + flag.
        ms = market_status.status(nm)
        mult = 1.0
        if ms:
            c["market_status"] = ms
            if ms["disqualifier"] == "hard":
                c["removed_reason"] = ms["label"]
                removed.append(c)
                continue
            mult *= ms["multiplier"]

        # 3. Mechanism crowding — me-too of the approved mechanism for this disease.
        mc = mechanism_crowding(nm, disease, occupied=occ, candidate_targets=tmap.get(nml))
        c["mechanism_crowding"] = mc
        mult *= mc["multiplier"]

        c["reality_multiplier"] = round(mult, 3)
        base = c.get("score") or c.get("composite_score") or 0
        try:
            c["discovery_rank_score"] = round(float(base) * mult, 4)
        except (TypeError, ValueError):
            c["discovery_rank_score"] = base
        leads.append(c)

    leads.sort(key=lambda x: x.get("discovery_rank_score") or 0, reverse=True)
    return {"leads": leads, "removed": removed}
