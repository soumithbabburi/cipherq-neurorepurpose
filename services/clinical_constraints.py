"""
Clinical Constraint Harmonizer (CCH)  —  the clinical-reality multiplier layer.
═══════════════════════════════════════════════════════════════════════════════
Multiplies a biological score by how well the drug fits the REAL-WORLD clinical
constraints of the disease:  R_final = R_bio × (S_tissue × S_duration × S_severity).

CONSOLIDATING, not a from-scratch rebuild:
  • S_tissue   — reuses CNS-MPO (barrier penetration). Narrow: only a HARD barrier
                 failure (CNS indication + clearly non-penetrant molecule) bites,
                 so it does not re-litigate the separately-shown developability axis.
  • drug risk  — reuses the FAERS serious-toxicity burden from services.safety_filter
                 (cached), so no new network and no re-implementation.
  • S_severity — NEW: a high-toxicity drug is only justified for a high-severity
                 disease; suggesting it for a benign/cosmetic condition is penalised.
  • S_duration — NEW: a high cumulative-toxicity drug is penalised for a CHRONIC
                 (lifelong) disease, but not for an ACUTE / life-threatening one.

Deliberately NOT double-counting the existing safety cross-filter (which penalises
organ-toxicity ↔ disease-ORGAN overlap): that is an organ-match check; CCH here is
a severity/chronicity-MATCH check — different axes.

FAIL-SOFT & HONEST: severity / chronicity are heuristic (keyword + therapeutic
area). Any factor we cannot determine stays 1.0 — we never invent a penalty from
missing data. Every penalty is bounded and flagged.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# ── Disease severity (proxy) ──────────────────────────────────────────────────
# High-severity / high-mortality contexts where a toxic drug is justified.
_HIGH_SEVERITY = (
    "cancer", "carcinoma", "neoplasm", "tumor", "tumour", "sarcoma", "leukemia",
    "leukaemia", "lymphoma", "myeloma", "melanoma", "glioma", "metasta",
    "failure", "fibrosis", "cirrhosis", "sepsis", "septic", "infarction",
    "ischemi", "ischaemi", "hemorrhage", "haemorrhage", "fatal", "terminal",
    "amyotrophic", "pulmonary hypertension", "cystic fibrosis",
    "amyloidosis", "muscular dystrophy", "aneurysm",
)
# Low-severity / benign / largely cosmetic — a high-toxicity drug is NOT justified.
_LOW_SEVERITY = (
    "vitiligo", "alopecia", "androgenetic", "hyperhidrosis", "rhinitis",
    "constipation", "hemorrhoid", "haemorrhoid", "acne", "rosacea",
    "seborrh", "wart", "dandruff", "halitosis", "cosmetic", "male pattern",
    "erectile", "dry eye", "presbyopia", "motion sickness", "jet lag",
)
_HIGH_SEVERITY_AREAS = ("cancer or benign tumor", "cardiovascular disorder")

# Hereditary cancer-PREDISPOSITION syndromes. Open Targets files these under an
# organ system (Cowden = "integumentary"), so the platform reads them as BENIGN —
# but they carry lifelong cancer risk and belong in the oncology-grade severity
# frame, not the benign-dermatology one. Recognized by named syndrome or by a
# predisposition/hereditary-cancer pattern.
_CANCER_PREDISPOSITION = (
    "cowden", "li-fraumeni", "li fraumeni", "lynch", "hamartoma tumor",
    "pten hamartoma", "peutz-jeghers", "peutz jeghers", "gorlin", "gardner",
    "familial adenomatous polyposis", "hereditary breast", "hereditary ovarian",
    "hereditary nonpolyposis", "von hippel", "von hippel-lindau", "birt-hogg",
    "carney complex", "multiple endocrine neoplasia", "retinoblastoma",
    "neurofibromatosis", "tuberous sclerosis", "xeroderma pigmentosum",
    "fanconi anemia", "ataxia telangiectasia", "bloom syndrome",
    "nijmegen breakage", "dyskeratosis congenita", "wilms",
    "cancer predisposition", "tumor predisposition", "tumour predisposition",
    "cancer susceptibility", "familial cancer",
)


def is_cancer_predisposition(disease_name: str) -> bool:
    """True for a hereditary cancer-predisposition / tumor syndrome."""
    return any(k in (disease_name or "").lower() for k in _CANCER_PREDISPOSITION)


def disease_severity(disease_name: str, areas: Optional[List[str]] = None) -> str:
    """'high' | 'low' | 'moderate' (moderate = unknown/neutral → no penalty)."""
    d = (disease_name or "").lower()
    blob = (d + " " + " ".join(areas or [])).lower()
    # A cancer-predisposition syndrome is oncology-grade severity, never benign,
    # regardless of the organ system Open Targets files it under.
    if is_cancer_predisposition(d):
        return "high"
    if any(k in d for k in _LOW_SEVERITY):
        return "low"
    if any(k in blob for k in _HIGH_SEVERITY) or any(a in blob for a in _HIGH_SEVERITY_AREAS):
        return "high"
    return "moderate"


# ── Disease chronicity (proxy) ────────────────────────────────────────────────
_ACUTE = ("acute", "infection", "sepsis", "septic", "poisoning", "overdose",
          "injury", "trauma", "hemorrhage", "haemorrhage", "wound", "burn",
          "emergency")
_CHRONIC = ("chronic", "syndrome", "disorder", "disease", "arthritis", "diabet",
            "sclerosis", "fibrosis", "hypertension", "degenerative")


# Localized / topically-treatable diseases — a systemic drug is a poor risk-benefit
# here when a LOCAL formulation could deliver the same target engagement without
# systemic toxicity (the vitiligo → oral-JAK vs ruxolitinib-cream problem).
_LOCALIZED = (
    "vitiligo", "alopecia", "psoriasis", "dermatit", "eczema", "acne", "rosacea",
    "skin", "cutaneous", "nail", "onychomycosis", "actinic keratosis",
    "conjunctivit", "dry eye", "keratit", "uveit", "glaucoma", "blepharit",
    "retinopath", "macular", "ocular", "otitis", "rhinitis", "nasal",
)


def disease_is_localized(disease_name: str) -> bool:
    d = (disease_name or "").lower()
    return any(k in d for k in _LOCALIZED)


_bw_cache: dict = {}


def has_boxed_warning(drug_name: str) -> bool:
    """True if the FDA label carries a BOXED (black-box) WARNING — the strongest
    systemic-risk signal, and one FAERS organ-toxicity misses (e.g. a JAK
    inhibitor's boxed warning for serious infections / MACE / malignancy). Cached,
    fail-soft → False."""
    d = (drug_name or "").strip().lower()
    if not d:
        return False
    if d in _bw_cache:
        return _bw_cache[d]
    out = False
    try:
        from services import http_client
        for field in ("generic_name", "brand_name", "substance_name"):
            r = http_client.get("https://api.fda.gov/drug/label.json", params={
                "search": f'openfda.{field}:"{d}"', "limit": 1}, timeout=10)
            if r is not None and r.ok:
                res = r.json().get("results", [])
                if res:
                    out = bool((res[0].get("boxed_warning") or [""])[0].strip())
                    break
    except Exception as e:
        logger.debug(f"boxed-warning lookup failed for {drug_name}: {e}")
    _bw_cache[d] = out
    return out


_ONCOLOGY_IND = ("cancer", "carcinoma", "neoplasm", "tumor", "tumour", "leukemia",
                 "leukaemia", "lymphoma", "myeloma", "melanoma", "sarcoma", "glioma")


def _is_oncology_drug(indications: str) -> bool:
    return any(k in (indications or "").lower() for k in _ONCOLOGY_IND)


_LOCAL_ROUTE_TERMS = ("topical", "ophthalmic", "otic", "cutaneous", "cream",
                      "ointment", "gel", "lotion", "transdermal", "nasal",
                      "inhalation", "intravitreal", "dermal")
_route_cache: dict = {}


def has_local_route(drug_name: str) -> bool:
    """True if the drug has an FDA-approved LOCAL formulation (topical, ophthalmic,
    inhaled…) — so it can treat a localized disease without systemic exposure. From
    the Orange Book DF;Route. Fail-soft → False."""
    d = (drug_name or "").strip().lower()
    if not d:
        return False
    if d in _route_cache:
        return _route_cache[d]
    found = False
    try:
        from services.orange_book import _rows, _ensure
        if _ensure():
            for row in _rows("products.txt"):
                ing = (row.get("Ingredient") or "").lower()
                if d in ing or ing.split()[0:1] == [d]:
                    route = (row.get("DF;Route") or "").lower()
                    if any(t in route for t in _LOCAL_ROUTE_TERMS):
                        found = True
                        break
    except Exception as e:
        logger.debug(f"route lookup failed for {drug_name}: {e}")
    _route_cache[d] = found
    return found


def disease_chronicity(disease_name: str) -> str:
    """'acute' | 'chronic' | 'unknown' (unknown → no duration penalty)."""
    d = (disease_name or "").lower()
    if any(k in d for k in _ACUTE):
        return "acute"
    if any(k in d for k in _CHRONIC):
        return "chronic"
    return "unknown"


# ── Drug risk (reuse FAERS serious-toxicity burden) ───────────────────────────
def drug_risk_profile(drug_name: str) -> Dict:
    """{risk 0..1, cytotoxic bool} from the FAERS serious-toxicity profile (reused
    from services.safety_filter, cached).

    `cytotoxic` = myelosuppression is a dominant serious signal — the hallmark of
    an antineoplastic / cytotoxic agent (palbociclib, doxorubicin). This is what
    distinguishes a drug that is genuinely disqualifying for a benign/chronic
    indication from one that merely has *a* notable organ AE it tolerates
    chronically (e.g. aspirin's GI bleeding). 0 / False when no data (no penalty)."""
    try:
        from services.safety_filter import drug_toxicity_profile
        tox = drug_toxicity_profile(drug_name) or {}
    except Exception as e:
        logger.debug(f"drug_risk failed for {drug_name}: {e}")
        tox = {}
    risk = round(min(1.0, sum(tox.values()) / 2.0), 3)
    cytotoxic = tox.get("hematologic", 0.0) >= 0.5
    return {"risk": risk, "cytotoxic": cytotoxic}


def drug_risk(drug_name: str) -> float:
    return drug_risk_profile(drug_name)["risk"]


# ── Tissue / barrier (reuse CNS-MPO, narrow hard-filter only) ─────────────────
_CNS_TERMS = ("brain", "cns", "central nervous", "neuro", "alzheimer", "parkinson",
              "multiple sclerosis", "epilep", "seizure", "dementia", "glioma",
              "encephal", "cerebral", "psychiat", "depression", "schizophren")


def _cns_barrier_factor(disease_name: str, smiles: str) -> tuple:
    """S_tissue for CNS indications only: a clearly non-penetrant molecule loses
    ground. Non-CNS or unknown MPO → 1.0."""
    d = (disease_name or "").lower()
    if not smiles or not any(t in d for t in _CNS_TERMS):
        return 1.0, None
    try:
        from services.cns_mpo import cns_mpo
        mpo = cns_mpo(smiles=smiles)
        score = mpo.get("score") if isinstance(mpo, dict) else None
        if score is None:
            return 1.0, None
        if score < 3.0:                     # clearly poor CNS exposure
            return 0.45, f"CNS indication but low brain penetration (CNS-MPO {score}/6)"
        if score < 4.0:
            return 0.8, None
        return 1.0, None
    except Exception:
        return 1.0, None


# ── Harmonizer ────────────────────────────────────────────────────────────────
# Softened 2026-07: raised from 0.1 so a single clinical-mismatch axis cannot crush
# an otherwise mechanistically-strong pair below the noise floor.
_FLOOR = 0.2


def harmonize(drug_name: str, disease_name: str, *,
              smiles: str = "", therapeutic_areas: Optional[List[str]] = None,
              has_trials: bool = False, indications: str = "") -> Dict:
    """R_bio multiplier from clinical constraints. Returns
    {multiplier, factors{tissue,duration,severity}, flags, penalized}. Fail-soft.

    has_trials: real clinical trials of this drug in this indication exist — humans
    have already judged the therapeutic window worth testing, so the toxicity
    penalties are lifted (empirical clinical movement overrides the heuristic)."""
    result = {"multiplier": 1.0, "penalized": False, "flags": [],
              "factors": {"tissue": 1.0, "duration": 1.0, "severity": 1.0}}
    if not drug_name or not disease_name:
        return result

    sev = disease_severity(disease_name, therapeutic_areas)
    chron = disease_chronicity(disease_name)
    rp = drug_risk_profile(drug_name)
    risk, cytotoxic = rp["risk"], rp["cytotoxic"]

    # A cytotoxic / antineoplastic agent is disqualifying outside high-severity
    # disease; a merely-notable-AE drug (risk based) gets a gentler ding.
    tox_word = "cytotoxic/myelosuppressive" if cytotoxic else "high-toxicity"

    # S_severity — the toxicity is only justified by disease severity. A cytotoxic
    # drug for a NON-high-severity disease is crushed; a high-risk (non-cytotoxic)
    # drug for a LOW-severity disease gets a proportional ding. HIGH severity → 1.0.
    s_sev = 1.0
    if sev != "high":
        if cytotoxic:
            s_sev = 0.30 if sev == "low" else 0.50          # softened (was 0.15/0.30)
        elif risk > 0.4 and sev == "low":
            s_sev = round(max(_FLOOR, 1.0 - 0.7 * risk), 3)
        if s_sev < 0.999:
            result["flags"].append(
                f"{tox_word.capitalize()} drug for a {sev}-severity indication - "
                f"therapeutic window unfavourable.")

    # S_duration — cumulative toxicity vs a CHRONIC (lifelong) regimen. Cytotoxic
    # agents are penalised hard for chronic use; acute disease → no penalty.
    s_dur = 1.0
    if chron == "chronic":
        if cytotoxic:
            s_dur = 0.50                                     # softened (was 0.30)
        elif risk > 0.5:
            s_dur = round(max(_FLOOR, 1.0 - 0.4 * risk), 3)
        if s_dur < 0.999:
            result["flags"].append(
                f"{tox_word.capitalize()} drug for a chronic (lifelong) indication "
                f"- long-term tolerability unlikely.")

    # S_route — localized vs systemic delivery. For a LOCALIZED, benign disease
    # (vitiligo, dermatitis, dry eye…) a systemic drug carrying real toxicity is a
    # poor risk-benefit when a LOCAL formulation could engage the same target
    # without systemic exposure. A drug WITH a topical/local formulation is spared
    # (its systemic penalties are waived — e.g. ruxolitinib CREAM for vitiligo);
    # an oral-only high-risk drug (upadacitinib, cerdulatinib) is penalised.
    s_route = 1.0
    if disease_is_localized(disease_name):
        if has_local_route(drug_name):
            # Local delivery avoids systemic toxicity — waive the systemic penalties.
            if s_sev < 0.9 or s_dur < 0.9:
                s_sev, s_dur = max(s_sev, 0.9), max(s_dur, 0.9)
                result["flags"] = [f for f in result["flags"]
                                   if "window" not in f and "tolerability" not in f]
                result["flags"].append(
                    "Local (topical/ophthalmic) formulation available - systemic "
                    "toxicity penalty waived for this localized indication.")
        elif cytotoxic or risk > 0.3 or has_boxed_warning(drug_name) or _is_oncology_drug(indications):
            s_route = 0.25
            result["flags"].append(
                "Systemic-only high-risk drug (boxed warning / oncology / cytotoxic) "
                "for a LOCALIZED indication - a local formulation would be required "
                "to justify the risk-benefit.")

    # S_tissue — CNS barrier hard-filter only.
    s_tis, tis_flag = _cns_barrier_factor(disease_name, smiles)
    if tis_flag:
        result["flags"].append(tis_flag)

    # Real trials in this indication override the toxicity heuristic — clinicians
    # already deemed the risk worth testing. Lift the toxicity floors (not the
    # tissue barrier, which is physical).
    if has_trials and (s_sev < 0.7 or s_dur < 0.7):
        s_sev, s_dur = max(s_sev, 0.7), max(s_dur, 0.7)
        result["flags"] = [f for f in result["flags"] if "therapeutic window" not in f
                           and "tolerability" not in f]
        result["flags"].append(
            "Toxicity penalty relaxed: real clinical trials of this drug in this "
            "indication already exist (human therapeutic-window judgment).")

    # Severity and duration are two views of ONE fact — "this drug is too toxic for
    # this clinical context" — so combine them with min() rather than multiplying,
    # which would double-count the same cytotoxicity (a cytotoxic drug for a chronic
    # moderate disease was getting 0.30×0.30≈0.09). Tissue (physical CNS barrier) and
    # route (local-vs-systemic delivery) are genuinely independent axes → still stack.
    _tox = min(s_sev, s_dur)
    mult = max(_FLOOR, round(_tox * s_tis * s_route, 3))
    result["factors"] = {"tissue": s_tis, "duration": s_dur, "severity": s_sev,
                         "route": s_route}
    result["multiplier"] = mult
    result["penalized"] = mult < 0.999
    result["disease_severity"] = sev
    result["disease_chronicity"] = chron
    result["drug_risk"] = risk
    result["cancer_predisposition"] = is_cancer_predisposition(disease_name)
    return result
