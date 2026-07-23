"""
Asset Acquisition Signal  (the flyer's 4-step engine, step 4)
═══════════════════════════════════════════════════════════════════════════════
Flags how readily the rights to a molecule can be OBTAINED for repurposing —
acquired, in-licensed, or freely developed via 505(b)(2) — and how attractive
the timing is. This is the BD layer the platform was missing: not "is this a
good scientific candidate" (the scorer answers that) but "can we get it, and
should we move now."

Every signal is grounded in authoritative FDA facts, not estimates:

  • Patent + regulatory-exclusivity expiry   → services.orange_book (active vs
    expired is a date fact). Drives "is it protected, and for how long."
  • Application types in the Orange Book      → NDA (Appl_Type 'N', originator/
    branded) vs ANDA (Appl_Type 'A', generic). ANDAs present ⇒ genericized ⇒
    freely sourceable. Count of ANDA sponsors ⇒ generic competition.
  • Marketing status (Type column)            → DISCN on every product ⇒ the
    originator has shelved the asset ⇒ classic in-licensing/acquisition target.
  • Reference Listed Drug (RLD = Yes)         → a 505(b)(2) bridge reference
    exists, so the de-risked regulatory pathway is open.

The output is an acquirability score in [0,1] (higher = easier + better-timed to
obtain for repurposing), a plain-English signal label, the regulatory pathway,
and the facts behind it.
"""
from __future__ import annotations

import csv
import logging
from datetime import date, datetime
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

try:                                              # package context (flask_app)
    from services.orange_book import orange_book_protection, OB_DIR
except ImportError:                               # script context
    from orange_book import orange_book_protection, OB_DIR  # type: ignore


def _parse_date(s: str) -> Optional[date]:
    s = (s or "").strip()
    for fmt in ("%b %d, %Y", "%Y-%m-%d", "%m/%d/%Y"):
        try:
            return datetime.strptime(s, fmt).date()
        except ValueError:
            continue
    return None


def _scan_products(name: str) -> Dict:
    """Originator, generic competition and marketing status from products.txt."""
    path = OB_DIR / "products.txt"
    out = {
        "found": False, "appl_types": set(), "nda_applicants": set(),
        "anda_applicants": set(), "statuses": set(), "rld": False,
        "earliest_approval": None, "originator": "",
    }
    if not path.exists():
        return out
    name_u = (name or "").strip().upper()
    if not name_u:
        return out
    with open(path, encoding="latin-1", newline="") as f:
        for row in csv.DictReader(f, delimiter="~"):
            ing = (row.get("Ingredient") or "").strip().upper()
            if not ing:
                continue
            # same matching policy as orange_book (exact or containment for >3 chars)
            tokens = [t.strip() for t in ing.split(";")]
            hit = any(name_u == t or (len(name_u) > 3 and (name_u in t or t in name_u))
                      for t in tokens)
            if not hit:
                continue
            out["found"] = True
            atype = (row.get("Appl_Type") or "").strip().upper()
            applicant = (row.get("Applicant_Full_Name") or row.get("Applicant") or "").strip()
            status = (row.get("Type") or "").strip().upper()
            out["appl_types"].add(atype)
            out["statuses"].add(status)
            if atype == "N" and applicant:
                out["nda_applicants"].add(applicant)
            elif atype == "A" and applicant:
                out["anda_applicants"].add(applicant)
            if (row.get("RLD") or "").strip().lower() == "yes":
                out["rld"] = True
            appd = _parse_date(row.get("Approval_Date") or "")
            if appd and (out["earliest_approval"] is None or appd < out["earliest_approval"]):
                out["earliest_approval"] = appd
                if atype == "N" and applicant:
                    out["originator"] = applicant
    if not out["originator"] and out["nda_applicants"]:
        out["originator"] = sorted(out["nda_applicants"])[0]
    return out


def _latest_active_expiry(protection: Dict) -> Optional[date]:
    latest = None
    for bucket in ("patents", "exclusivities"):
        for item in protection.get(bucket, []):
            if item.get("status") == "active" and item.get("expiry_iso"):
                d = _parse_date(item["expiry_iso"])
                if d and (latest is None or d > latest):
                    latest = d
    return latest


def _active_substance_expiry(protection: Dict) -> Optional[date]:
    """Latest ACTIVE *drug-substance* patent expiry — the only patent class that
    protects the molecule itself. Drug-product / method-of-use patents protect a
    specific formulation or indication and are precisely what a 505(b)(2) program
    is designed to engineer around, so they are NOT molecule-level barriers."""
    latest = None
    for p in protection.get("patents", []):
        if p.get("status") == "active" and "drug substance" in (p.get("type") or "") \
                and p.get("expiry_iso"):
            d = _parse_date(p["expiry_iso"])
            if d and (latest is None or d > latest):
                latest = d
    return latest


def _has_active_formulation_patents(protection: Dict) -> bool:
    """Any active patent that is product/method-of-use only (not drug substance)."""
    for p in protection.get("patents", []):
        t = (p.get("type") or "")
        if p.get("status") == "active" and "drug substance" not in t and \
                ("drug product" in t or "method-of-use" in t):
            return True
    return False


def acquisition_signal(drug_name: str) -> Dict:
    """Acquisition / in-licensing signal for a molecule, from FDA facts."""
    protection = orange_book_protection(drug_name)
    prod = _scan_products(drug_name)
    today = date.today()

    if not protection.get("available"):
        return {"available": False, "drug": drug_name,
                "signal": "Orange Book unavailable", "score": None}
    if not prod["found"] and not protection.get("patents") and not protection.get("exclusivities"):
        return {
            "available": True, "drug": drug_name, "score": 0.5,
            "signal": "Not in FDA Orange Book",
            "rationale": "No US NDA/ANDA record — likely a non-US, biologic, or "
                         "unapproved molecule. Acquirability cannot be inferred "
                         "from FDA data; assess via local registries.",
            "regulatory_pathway": "Unknown (no RLD)", "facts": {},
        }

    latest_expiry = _latest_active_expiry(protection)          # any active patent/excl
    subst_expiry = _active_substance_expiry(protection)        # molecule-level only
    has_formulation_patents = _has_active_formulation_patents(protection)
    molecule_protected = subst_expiry is not None
    subst_years = round((subst_expiry - today).days / 365.25, 1) if subst_expiry else 0.0

    has_generics = "A" in prod["appl_types"]
    n_competitors = len(prod["anda_applicants"])
    all_discontinued = bool(prod["statuses"]) and prod["statuses"] <= {"DISCN"}
    rld = prod["rld"]

    # The accessibility question is about the MOLECULE, not a branded line-extension.
    # Generic (ANDA) products on the market are proof the substance is freely
    # accessible — any later-dated patents protect specific formulations/uses,
    # which is exactly the white space a 505(b)(2) program occupies.
    formulation_note = (
        f" Residual formulation/method-of-use patents (to ~{latest_expiry.isoformat()}) "
        "protect specific branded products, not the molecule — a 505(b)(2) program "
        "is designed around them."
        if has_formulation_patents and latest_expiry else "")
    crowd_note = (f" Note: {n_competitors} generic sponsors — a commoditised market, so "
                  "differentiation (new indication/formulation) is essential."
                  if n_competitors >= 15 else "")

    # ── Classification (ordered most → least actionable) ──────────────────────
    if all_discontinued:
        signal, score = "Shelved / discontinued — in-licensing target", 0.92
        why = ("Every Orange Book product for this molecule is marked "
               "discontinued (DISCN); the originator has stopped marketing it. "
               "Rights are typically available to in-license or acquire — the "
               "canonical repurposing opportunity.")
    elif has_generics:
        # Molecule is genericized ⇒ freely developable regardless of line-extension patents.
        if has_formulation_patents:
            signal, score = "Genericized base — 505(b)(2) white space", 0.80
            why = (f"{n_competitors or 'Multiple'} generic (ANDA) sponsor(s) confirm the "
                   "drug substance is off-patent and freely sourceable." + formulation_note + crowd_note)
        else:
            signal, score = "Fully off-patent — open for 505(b)(2)", 0.85
            why = (f"{n_competitors or 'Multiple'} generic sponsor(s) and no active "
                   "substance patent — develop a differentiated product directly, "
                   "no license required." + crowd_note)
    elif molecule_protected and subst_years <= 3:
        signal, score = f"Substance patent expires in ~{subst_years}y — in-licensing window", 0.70
        why = (f"The drug-substance patent lapses ~{subst_expiry.isoformat()}. Prime "
               "window to in-license now and launch a repurposed product as the "
               "molecule opens to generics.")
    elif molecule_protected and subst_years <= 6:
        signal, score = f"Originator-protected molecule (~{subst_years}y) — partnership", 0.50
        why = (f"Drug-substance patent active to ~{subst_expiry.isoformat()}. A "
               "new-indication program is viable but needs a partnership/license with "
               f"the originator ({prod['originator'] or 'unknown'}).")
    elif molecule_protected:
        signal, score = f"Originator-protected molecule ({subst_years}y runway) — acquisition only", 0.35
        why = (f"Long drug-substance protection (to ~{subst_expiry.isoformat()}). Access "
               "requires acquisition or a licensing deal with the originator "
               f"({prod['originator'] or 'unknown'}).")
    elif has_formulation_patents:
        # No generics yet, no substance patent — molecule open, only formulation patents.
        signal, score = "No generics; only formulation patents — 505(b)(2) white space", 0.82
        why = ("The drug substance is unprotected and no generic has entered yet — "
               "first-mover white space." + formulation_note)
    else:
        signal, score = "Off-patent, no generics — first-mover white space", 0.85
        why = ("All protection has lapsed and no generic has entered. Freely "
               "developable with first-mover advantage on a new indication.")

    # Slight uplift when a 505(b)(2) bridge reference exists (de-risked pathway).
    if rld and score is not None and score < 0.9:
        score = round(min(0.95, score + 0.05), 2)

    pathway = ("505(b)(2) — bridge to RLD" if rld
               else "505(b)(2) likely (no RLD flagged)" if (has_generics or not molecule_protected)
               else "Full NDA or license-dependent")

    return {
        "available": True, "drug": drug_name,
        "signal": signal, "score": score,
        "rationale": why,
        "regulatory_pathway": pathway,
        "facts": {
            "molecule_protected": molecule_protected,
            "substance_patent_until": subst_expiry.isoformat() if subst_expiry else None,
            "substance_years_to_open": subst_years,
            "latest_any_protection_until": latest_expiry.isoformat() if latest_expiry else None,
            "has_active_formulation_patents": has_formulation_patents,
            "has_generics": has_generics,
            "generic_competitors": n_competitors,
            "all_discontinued": all_discontinued,
            "rld_available": rld,
            "originator": prod["originator"],
            "active_patents": sum(1 for p in protection.get("patents", []) if p.get("status") == "active"),
            "active_exclusivities": sum(1 for e in protection.get("exclusivities", []) if e.get("status") == "active"),
        },
        "source": "FDA Orange Book (patents, exclusivities, application types, marketing status)",
    }


def rank_for_acquisition(drug_names: List[str]) -> List[Dict]:
    """Score a list of molecules, best acquisition signal first."""
    out = [acquisition_signal(d) for d in drug_names]
    return sorted(out, key=lambda r: (r.get("score") or 0), reverse=True)


if __name__ == "__main__":
    for d in ["thalidomide", "metformin", "sildenafil", "apixaban", "minoxidil"]:
        r = acquisition_signal(d)
        f = r.get("facts", {})
        print(f"\n{d.upper():<14} score={r.get('score')}  {r.get('signal')}")
        print(f"   pathway: {r.get('regulatory_pathway')}")
        if f:
            print(f"   substance_patent={f.get('substance_patent_until')} "
                  f"formulation_patents={f.get('has_active_formulation_patents')} "
                  f"generics={f.get('has_generics')}({f.get('generic_competitors')}) "
                  f"discontinued={f.get('all_discontinued')} rld={f.get('rld_available')}")
