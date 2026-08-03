"""
505(b)(2) Feasibility & Exclusivity (Layer 4).
════════════════════════════════════════════════════════════════════════════════
505(b)(2) is a regulatory/business STRATEGY, but a target-matching score is blind to
it. This layer turns "is this a fundable repurposing filing?" into a real calculation:

  • EXCLUSIVITY RUNWAY — how long the BASE compound is still protected (FDA Orange Book
    patents/exclusivity for small molecules; FDA Purple Book biosimilar cliff for
    biologics), split into PRIMARY (drug-substance) vs SECONDARY (formulation / method-
    of-use) patents.
  • REGULATORY ROUTE — 505(b)(2) for a small molecule; an sBLA / label expansion for a
    biologic (505(b)(2) is small-molecule only; biosimilar entry is via 351(k)).
  • FORMULATION OPPORTUNITY — the highest-value angle: reformulating for a compartment-
    alized indication (an IV/systemic biologic → inhaled / subcutaneous / topical) mints
    NEW formulation + device IP and a differentiated filing, turning a low-exclusivity
    asset into a fundable one. Reuses services.delivery_feasibility.
  • ORPHAN — 7-year orphan exclusivity strengthens a small-market case.

Grounded in real FDA data; fail-soft (unknown, not fabricated, when data is absent).
"""
from __future__ import annotations

import logging
from datetime import date
from typing import Dict, Optional

logger = logging.getLogger(__name__)


def _years_until(iso: Optional[str]) -> Optional[float]:
    try:
        return round((date.fromisoformat(iso) - date.today()).days / 365.25, 1)
    except Exception:
        return None


def exclusivity_profile(drug_name: str, is_biologic: bool = False) -> Dict:
    """Drug-level base-compound protection (constant across a drug's candidate indications).
    Compute ONCE per drug."""
    out = {"runway_years": None, "route": None, "patents_primary": 0,
           "patents_secondary": 0, "latest_active": None, "source": None, "narrative": ""}

    if is_biologic:
        out["route"] = ("sBLA / label expansion (biologic — 505(b)(2) is small-molecule "
                        "only; biosimilar entry is via the 351(k) pathway)")
        try:
            from services.purple_book import biologic_status
            b = biologic_status(drug_name) or {}
            if b.get("found"):
                out["source"] = "FDA Purple Book"
                cliff = b.get("cliff")
                out["latest_active"] = cliff
                out["runway_years"] = _years_until(cliff) if cliff else None
                out["n_biosimilars"] = b.get("n_biosimilars", 0)
                if b.get("cliff_open"):
                    out["narrative"] = (f"Biosimilar cliff is open ({b.get('n_biosimilars', 0)} "
                                        f"biosimilars) — base exclusivity has lapsed, so a new "
                                        f"indication needs its own method-of-use / formulation IP.")
                elif cliff:
                    out["narrative"] = (f"Base biologic protected until ~{cliff} "
                                        f"({out['runway_years']}y); a new indication extends the "
                                        f"franchise under the existing BLA.")
        except Exception as e:
            logger.debug(f"purple book lookup: {e}")
        return out

    out["route"] = ("505(b)(2) — references the approved drug's safety/PK; a new-indication "
                    "study earns 3-year exclusivity")
    # Two DISTINCT protections, never conflated (audited 2026-07): FDA regulatory
    # exclusivity (NCE 5y / orphan 7y / new-indication 3y — a statutory market block that
    # bars ANDA/505(b)(2) reliance) vs PATENT life (which may run a decade-plus longer on a
    # late formulation/method-of-use patent). The headline "exclusivity runway" reports the
    # REGULATORY exclusivity; patent life is reported separately so a 2044 formulation patent
    # is never presented as "17.9y of exclusivity". The FRANCHISE runway (the longer of the
    # two — competitors are blocked while EITHER holds) drives the feasibility SCORE only.
    out.update({"patent_expiry": None, "patent_runway_years": None,
                "regulatory_exclusivity_expiry": None, "regulatory_exclusivity_code": None,
                "regulatory_exclusivity_label": None, "regulatory_runway_years": None,
                "franchise_runway_years": None})
    try:
        from services.orange_book import orange_book_protection
        ob = orange_book_protection(drug_name) or {}
        if ob.get("available") and (ob.get("patents") or ob.get("exclusivities")):
            out["source"] = "FDA Orange Book"
            active = [p for p in ob.get("patents", []) if p.get("status") == "active"]
            out["patents_primary"] = sum(1 for p in active if "drug substance" in (p.get("type") or ""))
            out["patents_secondary"] = sum(1 for p in active if "drug substance" not in (p.get("type") or ""))
            pat_exps = [p["expiry_iso"] for p in active if p.get("expiry_iso")]
            latest_patent = max(pat_exps) if pat_exps else None
            out["patent_expiry"] = latest_patent
            out["patent_runway_years"] = _years_until(latest_patent) if latest_patent else None

            # Regulatory exclusivity — the statutory market block (NCE/ODE/NP/I…), distinct
            # from patent life. Take the latest-expiring ACTIVE exclusivity as the headline.
            active_excl = [e for e in ob.get("exclusivities", [])
                           if e.get("status") == "active" and e.get("expiry_iso")]
            reg = max(active_excl, key=lambda e: e["expiry_iso"]) if active_excl else None
            if reg:
                out["regulatory_exclusivity_expiry"] = reg["expiry_iso"]
                out["regulatory_exclusivity_code"] = reg.get("code")
                out["regulatory_exclusivity_label"] = reg.get("label")
                out["regulatory_runway_years"] = _years_until(reg["expiry_iso"])

            # Headline "exclusivity runway" = REGULATORY exclusivity (NOT the latest patent).
            out["runway_years"] = out["regulatory_runway_years"]
            out["latest_active"] = out["regulatory_exclusivity_expiry"]
            # Franchise runway (feasibility scoring only) = the longer of patent vs exclusivity.
            _runways = [r for r in (out["patent_runway_years"], out["regulatory_runway_years"])
                        if r is not None]
            out["franchise_runway_years"] = max(_runways) if _runways else None

            _reg_bit = (f"Regulatory exclusivity ({out['regulatory_exclusivity_label']}) to ~"
                        f"{out['regulatory_exclusivity_expiry']} ({out['regulatory_runway_years']}y)"
                        if reg else "No unexpired FDA regulatory exclusivity")
            _pat_bit = (f"latest patent to ~{latest_patent} ({out['patent_runway_years']}y; "
                        f"{out['patents_primary']} substance / {out['patents_secondary']} "
                        f"formulation-or-use)" if latest_patent else "no active Orange Book patents")
            out["narrative"] = (f"{_reg_bit}; {_pat_bit}. Exclusivity and patent life are distinct — "
                                f"the franchise is blocked while either holds; repurposing can add its "
                                f"own new-indication (3y) exclusivity + method-of-use IP.")
            if not latest_patent and not reg:
                out["narrative"] = ("Base compound appears off-patent with no active exclusivity — a "
                                    "repurposed use needs NEW method-of-use / formulation IP; a "
                                    "505(b)(2) reformulation is the play.")
        else:
            out["narrative"] = ("No Orange Book patents or exclusivity on file — likely off-patent; "
                                "feasibility rests on new-indication (3y) exclusivity + any formulation IP.")
    except Exception as e:
        logger.debug(f"orange book lookup: {e}")
    return out


def feasibility(disease: str, drug_name: str, excl: Dict, is_biologic: bool = False,
                smiles: str = "", disease_value: Optional[Dict] = None) -> Dict:
    """Per-candidate 505(b)(2) feasibility: base runway + regulatory route + FORMULATION
    opportunity (disease-specific) + orphan. Returns a 20–95 score + components + flags."""
    flags = []
    score = 55                                   # a 505(b)(2)/sBLA path always exists for an approved drug
    excl = excl or {}

    # Formulation opportunity — the high-value angle (localized reformulation → new IP).
    form = None
    try:
        from services.delivery_feasibility import assess_delivery
        d = assess_delivery(disease, smiles=smiles, drug_name=drug_name) or {}
        route = d.get("route") or ""
        if d.get("localized") and d.get("deliverable") and route and "systemic" not in route.lower():
            form = {"opportunity": True, "route": route, "compartment": d.get("compartment"),
                    "note": (f"Reformulation opportunity — a {route} formulation for this "
                             f"compartmentalized indication mints NEW formulation/device IP and a "
                             f"differentiated 505(b)(2) filing (vs the systemic base).")}
            score += 22
            flags.append(form["note"])
    except Exception as e:
        logger.debug(f"delivery feasibility: {e}")

    runway = excl.get("runway_years")                      # headline = REGULATORY exclusivity
    # Feasibility SCORE uses the FRANCHISE runway (longer of patent vs exclusivity) — a drug
    # blocked by a long patent is still a defensible franchise even after regulatory
    # exclusivity lapses. Falls back to the regulatory runway when patent life is unknown.
    franchise = excl.get("franchise_runway_years")
    if franchise is None:
        franchise = runway
    if franchise is not None and franchise >= 5:
        score += 12
        flags.append(f"Base franchise protected ~{franchise}y (patent or exclusivity) — "
                     f"repurposing extends it, not a generic race.")
    elif franchise is not None and franchise < 1:
        score -= 6
        flags.append("Base compound is generic — exclusivity must come from the new use / formulation.")

    if (disease_value or {}).get("is_orphan"):
        score += 10
        flags.append("Orphan-eligible — 7-year orphan exclusivity strengthens the case.")
    if is_biologic:
        score += 6                               # own-BLA label expansion is a clean regulatory path

    score = int(max(20, min(95, score)))
    return {"feasibility_score": score, "route": excl.get("route"), "runway_years": runway,
            "patents_primary": excl.get("patents_primary"), "patents_secondary": excl.get("patents_secondary"),
            # Distinct patent vs regulatory-exclusivity fields (never conflated in display).
            "patent_expiry": excl.get("patent_expiry"),
            "patent_runway_years": excl.get("patent_runway_years"),
            "regulatory_exclusivity_expiry": excl.get("regulatory_exclusivity_expiry"),
            "regulatory_exclusivity_label": excl.get("regulatory_exclusivity_label"),
            "regulatory_runway_years": excl.get("regulatory_runway_years"),
            "franchise_runway_years": excl.get("franchise_runway_years"),
            "formulation_opportunity": form, "exclusivity_source": excl.get("source"),
            "exclusivity_narrative": excl.get("narrative"), "flags": flags,
            "tier": ("Fundable" if score >= 75 else "Workable" if score >= 55 else "Thin")}
