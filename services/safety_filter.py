"""
Safety / toxicity cross-filter  —  generalized negative-safety scoring.
═══════════════════════════════════════════════════════════════════════════════
A repurposing match is only credible if the drug is not actively dangerous to the
ORGAN the disease lives in. This module overlays a drug's real-world serious
adverse-event profile (FDA FAERS via openFDA) onto the target disease's organ
pathology and applies a NEGATIVE multiplier when they collide.

The logic is drug/disease-AGNOSTIC — it reasons over a shared organ-system
taxonomy, not hardcoded pairs:

    disease → implicated organ system(s)  (+ "critical" if the disease IS
                                           organ failure / damage / fibrosis)
    drug    → organ system(s) it is toxic to  (from SERIOUS FAERS reactions)
    collide on a critical organ            → strong penalty  (×~0.15)
    collide on an involved organ           → mild penalty    (×~0.55)

Two-tier taxonomy — deliberately asymmetric:
  • DISEASE side is BROAD: any pulmonary disease (asthma, COPD, fibrosis) maps to
    the "pulmonary" organ.
  • DRUG-TOXICITY side is STRICT: only genuine SERIOUS organ-damage MedDRA terms
    count (pulmonary embolism, pneumonitis, hepatic failure…), NOT nonspecific
    symptoms (dyspnoea, nausea, fatigue) that would falsely implicate an organ.

Canonical example this fixes: an anti-angiogenic that carries a serious
pulmonary-embolism / thrombosis signal is down-ranked for a pulmonary indication
(e.g. idiopathic pulmonary fibrosis), automatically and generally.

FAIL-SOFT & HONEST: if no serious adverse-event data can be resolved for the
drug, the multiplier is exactly 1.0 — we never fabricate a penalty from missing
data, and never block scoring on a network hiccup.
"""
from __future__ import annotations

import re
import logging
from typing import Dict, List, Tuple

from services import http_client

logger = logging.getLogger(__name__)

OPENFDA_EVENT = "https://api.fda.gov/drug/event.json"

# ── DISEASE → organ taxonomy (BROAD) ──────────────────────────────────────────
# Substring stems matched against a disease NAME. Broad on purpose: a disease
# named for an organ belongs to that organ system.
_DISEASE_ORGAN_STEMS: Dict[str, Tuple[str, ...]] = {
    "pulmonary":        ("pulmonary", "lung", "respiratory", "pneumon", "bronch",
                         "pleural", "alveol", "interstitial lung", "copd",
                         "asthma", "ards", "cystic fibrosis"),
    "cardiac":          ("cardiac", "myocard", "cardiomyopath", "heart",
                         "arrhythmi", "atrial fibrill", "coronary", "angina",
                         "pericard", "valv", "ejection fraction"),
    "vascular":         ("vascular", "thrombos", "thromboembol", "embolism",
                         "atheroscler", "aneurysm", "ischemic", "ischaemic",
                         "stroke", "peripheral arter", "venous", "vasculit",
                         "hypertension", "raynaud"),
    "hepatic":          ("hepat", "liver", "cirrhos", "cholestas", "biliary"),
    "renal":            ("renal", "kidney", "nephro", "nephrit", "glomerul"),
    "hematologic":      ("anemia", "anaemia", "neutropen", "thrombocytopen",
                         "leukemi", "leukaemi", "lymphom", "myelodysplas",
                         "myelofibros", "hematolog", "haematolog", "aplast"),
    "cns":              ("epilep", "seizure", "encephal", "cerebral",
                         "cognit", "dementia", "parkinson", "multiple sclerosis",
                         "neurodegener", "meningit", "alzheimer"),
    "gastrointestinal": ("colitis", "pancreatit", "gastrointestinal", "bowel",
                         "crohn", "enterocolitis", "gastric", "intestinal",
                         "esophag", "oesophag", "gastritis", "ulcerative",
                         "coeliac", "celiac", "diverticul", "eosinophilic gastro"),
    "dermatologic":     ("psorias", "dermat", "cutaneous", "epidermal",
                         "skin", "vitiligo", "eczema", "urticaria", "acne",
                         "pemphig", "ichthyos"),
    "ocular":           ("retinopath", "macular", "ocular", "retinal", "glaucoma",
                         "uveitis", "optic", "keratit", "conjunctiv"),
    "endocrine":        ("thyroid", "adrenal", "endocrin", "pituitary", "diabet",
                         "metabolic syndrome", "obesity", "dyslipid", "hyperlipid",
                         "hypercholester", "acromegaly", "cushing"),
    # Musculoskeletal / connective tissue — previously unmapped, so rheumatoid
    # arthritis, myopathies, and myositis fell through with no organ and the
    # organ-toxicity check silently skipped (a myotoxic drug went unpenalised).
    "musculoskeletal":  ("arthritis", "arthropath", "arthr", "synovial", "spondyl",
                         "myopath", "myositis", "myasthen", "muscular dystrophy",
                         "rhabdomyol", "osteoarthr", "osteoporos", "gout",
                         "tendin", "fibromyalg", "polymyalg", "sarcopeni",
                         "rheumatoid", "connective tissue", "lupus", "scleroderma",
                         "sjogren", "sjögren", "sarcoid", "polymyositis",
                         "dermatomyositis", "still's disease", "still disease"),
}

# ── DRUG-TOXICITY → organ taxonomy (STRICT) ───────────────────────────────────
# Substring stems matched against SERIOUS FAERS MedDRA reaction terms. Only
# genuine organ DAMAGE / failure / structural toxicity — NOT nonspecific symptoms
# (dyspnoea, cough, nausea, fatigue, pain), which are excluded by omission.
_TOX_ORGAN_STEMS: Dict[str, Tuple[str, ...]] = {
    "pulmonary":        ("pulmonary embolism", "pulmonary fibros",
                         "pulmonary hypertension", "pulmonary haemorrhage",
                         "pulmonary oedema", "pulmonary edema", "pneumonitis",
                         "interstitial lung", "respiratory failure",
                         "acute respiratory distress", "lung disorder",
                         "pleural effusion"),
    "cardiac":          ("cardiac failure", "myocardial infarct", "cardiomyopath",
                         "cardiac arrest", "arrhythmi", "atrial fibrill",
                         "ventricular", "qt prolong", "torsade", "myocarditis",
                         "cardiotox", "pericardi", "ejection fraction"),
    "vascular":         ("thrombos", "thromboembol", "embolism", "infarct",
                         "ischaemi", "ischemi", "cerebrovascular", "haemorrhage",
                         "hemorrhage", "aneurysm", "vasculit", "hypertensive crisis"),
    "hepatic":          ("hepatotox", "hepatic failure", "hepatitis",
                         "hepatic necros", "liver injury", "hepatic function",
                         "hyperbilirubin", "cholestas", "transaminases increased"),
    "renal":            ("renal failure", "renal impairment", "acute kidney",
                         "nephrotox", "nephritis", "renal insufficien",
                         "renal necros"),
    "hematologic":      ("neutropen", "thrombocytopen", "agranulocyt",
                         "pancytopen", "aplastic anaemia", "aplastic anemia",
                         "bone marrow", "leukopen", "febrile neutropen",
                         "myelosuppress", "haemolytic", "hemolytic"),
    "cns":              ("seizure", "convuls", "encephalopath", "intracranial",
                         "cerebral haemorrhage", "cerebral hemorrhage",
                         "posterior reversible", "stroke"),
    "gastrointestinal": ("gastrointestinal haemorrhage", "gastrointestinal hemorrhage",
                         "gastrointestinal perforation", "colitis", "pancreatitis",
                         "ileus", "intestinal perforation", "gastric perforation"),
    "dermatologic":     ("stevens-johnson", "toxic epidermal", "epidermal necrolysis",
                         "erythema multiforme", "drug reaction with eosinophilia"),
    "ocular":           ("retinopath", "retinal", "optic neuritis", "vision loss",
                         "blindness", "retinal vein"),
    "endocrine":        ("adrenal insufficien", "thyroiditis", "hypophysitis",
                         "hypothyroidism"),
    # Genuine structural muscle toxicity only (not nonspecific "muscle pain"/weakness).
    "musculoskeletal":  ("rhabdomyolysis", "myopathy", "myositis", "myopathy toxic",
                         "necrotising myopathy", "necrotizing myopathy",
                         "muscle necrosis", "osteonecrosis"),
}

# Tokens in a DISEASE name that mean the disease IS organ failure / structural
# damage — a toxicity to that organ is then especially unacceptable ("critical").
_CRITICAL_DISEASE_MARKERS = (
    "failure", "fibros", "insufficien", "cirrhos", "infarct", "injury",
    "end-stage", "end stage", "damage", "necros", "ischemi", "ischaemi",
    "arrest", "acute respiratory", "decompensat",
)


# ── Disease → organ profile ───────────────────────────────────────────────────
def disease_organ_profile(disease_name: str) -> Dict[str, str]:
    """{organ_system: 'critical' | 'involved'} implicated by a disease name."""
    d = (disease_name or "").lower()
    if not d:
        return {}
    is_critical = any(m in d for m in _CRITICAL_DISEASE_MARKERS)
    out: Dict[str, str] = {}
    for organ, stems in _DISEASE_ORGAN_STEMS.items():
        if any(s in d for s in stems):
            out[organ] = "critical" if is_critical else "involved"
    return out


def therapeutic_organs(indications_text: str) -> set:
    """Organ systems a drug is THERAPEUTICALLY developed against, read from its
    indications text. Layer 3 (Safety Contextualization): a serious-AE signal in one of
    THESE organs is confounded by the therapeutic context — the drug is MEANT to act
    there, and its FAERS reports in that organ come from the treated disease population,
    not from off-target damage. So a respiratory drug should NOT be safety-penalised for a
    respiratory disease (mepolizumab→ABPA). Distinct from true off-therapeutic toxicity
    (a hepatotoxic drug NOT developed for any liver disease, used in a liver disease)."""
    organs = set()
    for chunk in re.split(r"[;,/|]| and ", indications_text or ""):
        organs |= set(disease_organ_profile(chunk).keys())
    return organs


# ── Drug → toxicity organ profile (SERIOUS FAERS reactions) ───────────────────
def _classify_tox_term(term: str) -> List[str]:
    t = (term or "").lower()
    return [organ for organ, stems in _TOX_ORGAN_STEMS.items()
            if any(s in t for s in stems)]


_syn_cache: dict = {}


def _name_variants(drug: str) -> List[str]:
    """[drug] + its trade names, so a FAERS lookup that misses on the generic name
    still finds reports filed under the BRAND (e.g. palbociclib → IBRANCE). FAERS
    keys many reports by medicinalproduct=brand, so the generic alone silently
    returns nothing — the source of a phantom risk=0 for whole drug classes."""
    d = (drug or "").strip()
    if not d:
        return []
    key = d.lower()
    if key in _syn_cache:
        return _syn_cache[key]
    variants = [d]
    # Brand names from the FDA label (openFDA maps generic↔brand reliably; ChEMBL
    # synonyms give salt/research-code forms, not the trade name FAERS keys on).
    try:
        for field in ("generic_name", "substance_name"):
            r = http_client.get("https://api.fda.gov/drug/label.json", params={
                "search": f'openfda.{field}:"{d}"', "limit": 1}, timeout=10)
            if r is not None and r.ok:
                res = r.json().get("results", [])
                if res:
                    for bn in (res[0].get("openfda", {}) or {}).get("brand_name", []) or []:
                        if bn and bn.lower() not in {v.lower() for v in variants}:
                            variants.append(bn)
                    break
    except Exception as e:
        logger.debug(f"brand-name fetch failed for {drug}: {e}")
    variants = variants[:6]                       # cap the FAERS retries
    _syn_cache[key] = variants
    return variants


def _faers_one(name: str, n: int) -> List[Tuple[str, int]]:
    """Serious FAERS reactions for ONE exact name string.
    NOTE: the ` AND serious:1` filter must use spaces (encoded to +), NOT a
    literal `+AND+` — the latter 404s through the request encoder."""
    d = name.strip().lower()
    search = (f'(patient.drug.openfda.generic_name:"{d}"'
              f'+patient.drug.openfda.brand_name:"{d}"'
              f'+patient.drug.medicinalproduct:"{d}") AND serious:1')
    try:
        r = http_client.get(OPENFDA_EVENT, params={
            "search": search,
            "count": "patient.reaction.reactionmeddrapt.exact"}, timeout=12)
        if r is not None and r.ok:
            return [(x["term"], int(x["count"]))
                    for x in r.json().get("results", [])[:n] if x.get("term")]
    except Exception as e:
        logger.debug(f"FAERS fetch failed for {name}: {e}")
    return []


def _faers_serious_reactions(drug: str, n: int = 50) -> List[Tuple[str, int]]:
    """Top SERIOUS FAERS reaction terms for a drug — tries the generic name, then
    its brand synonyms, so a brand-keyed profile is not missed. [] on total miss."""
    if not (drug or "").strip():
        return []
    for name in _name_variants(drug):
        rx = _faers_one(name, n)
        if rx:
            return rx
    return []


def _organ_shares(reactions: List[Tuple[str, int]]) -> Dict[str, float]:
    """{organ: share of serious-AE mass} via the STRICT organ-damage taxonomy."""
    total = sum(c for _, c in reactions) or 1
    mass: Dict[str, int] = {}
    for term, cnt in reactions:
        for organ in _classify_tox_term(term):
            mass[organ] = mass.get(organ, 0) + cnt
    return {o: m / total for o, m in mass.items()}


# Disproportionality thresholds. A drug's organ toxicity only counts to the extent
# it is OVER-REPRESENTED vs the all-drug background — this is the standard
# pharmacovigilance correction (a proportional-reporting-ratio style lift) and it
# is what removes confounding-by-indication (e.g. aspirin's serious reports are
# full of cardiac/vascular events because cardiac patients take aspirin — that
# signal sits at ~1.8× background and is correctly ignored, while aspirin's GI
# bleeding at ~5× background is kept).
_LIFT_MIN, _LIFT_FULL = 2.0, 5.0        # lift<2 → ignored; ≥5 → full strength
_BG_FLOOR = 0.002                       # min background share (avoid div blow-up)
_MIN_DRUG_SHARE = 0.01                  # organ must be ≥1% of the drug's serious mass

_bg_cache: Dict[str, float] = {}
_bg_loaded = False


def _background_organ_shares() -> Dict[str, float]:
    """All-drug organ shares among SERIOUS reports — the disproportionality
    denominator. Fetched once, cached; {} on failure (→ lift disabled, fail-open
    to raw shares is NOT done — no background means no penalty)."""
    global _bg_loaded, _bg_cache
    if _bg_loaded:
        return _bg_cache
    _bg_loaded = True
    try:
        r = http_client.get(OPENFDA_EVENT, params={
            "search": "serious:1",
            "count": "patient.reaction.reactionmeddrapt.exact"}, timeout=15)
        if r is not None and r.ok:
            rows = [(x["term"], int(x["count"]))
                    for x in r.json().get("results", [])[:100] if x.get("term")]
            _bg_cache = _organ_shares(rows)
    except Exception as e:
        logger.debug(f"FAERS background fetch failed: {e}")
    return _bg_cache


def drug_toxicity_profile(drug: str) -> Dict[str, float]:
    """{organ_system: weight 0..1} of a drug's SERIOUS organ toxicity, corrected
    for background reporting (disproportionality). weight ramps 0→1 as the drug's
    organ share rises from 2× to 5× the all-drug background share. Empty dict when
    no data (→ no penalty)."""
    reactions = _faers_serious_reactions(drug)
    if not reactions:
        return {}
    drug_shares = _organ_shares(reactions)
    bg = _background_organ_shares()
    if not bg:                              # no denominator → no penalty (honest)
        return {}
    profile: Dict[str, float] = {}
    for organ, share in drug_shares.items():
        if share < _MIN_DRUG_SHARE:
            continue
        lift = share / max(bg.get(organ, 0.0), _BG_FLOOR)
        if lift <= _LIFT_MIN:
            continue
        w = min(1.0, (lift - _LIFT_MIN) / (_LIFT_FULL - _LIFT_MIN))
        if w > 0:
            profile[organ] = round(w, 3)
    return profile


# ── Cross-filter ──────────────────────────────────────────────────────────────
_CRITICAL_FLOOR = 0.15      # strongest allowed knock-down for a critical clash
_INVOLVED_FLOOR = 0.55      # milder knock-down for a non-critical clash
_MULTIPLIER_FLOOR = 0.10    # never zero a score out entirely


def assess(drug_name: str, disease_name: str, therapeutic_organs=None) -> Dict:
    """Cross a drug's serious-toxicity organ profile against a disease's organ
    pathology. Returns a bounded penalty multiplier, human-readable flags, and the
    underlying profiles for display. Fail-soft: multiplier 1.0 when data is absent.

    Layer 3 — Safety Contextualization: `therapeutic_organs` is the set of organ systems
    the drug is DEVELOPED to act on (from its indications). A toxicity signal in one of
    those organs is confounded by the therapeutic context (a lung drug's respiratory FAERS
    events come from its asthma population, not lung damage), so it is NOT penalised —
    only a toxicity in an organ the drug is NOT therapeutic for is a real contraindication.
    """
    result = {"multiplier": 1.0, "penalized": False, "flags": [],
              "drug_toxicity": {}, "disease_organs": {}}
    disease_organs = disease_organ_profile(disease_name)
    if not disease_organs:
        return result
    tox = drug_toxicity_profile(drug_name)
    if not tox:                                   # no data → no penalty (honest)
        result["disease_organs"] = disease_organs
        return result

    result["drug_toxicity"] = tox
    result["disease_organs"] = disease_organs
    therapeutic = set(therapeutic_organs or [])

    multiplier = 1.0
    for organ, severity in disease_organs.items():
        w = tox.get(organ, 0.0)
        if not w:
            continue
        if organ in therapeutic:
            # therapeutic-organ overlap — the drug is developed to act here; its FAERS
            # signal in this organ is the treated-disease context, not a contraindication.
            result.setdefault("suppressed_organs", []).append(organ)
            result["flags"].append(
                f"{organ.replace('_','/')} signal is therapeutic-organ overlap "
                f"(drug is developed to act on this organ) — not a contraindication, no penalty.")
            continue
        if severity == "critical":
            factor = 1.0 - (1.0 - _CRITICAL_FLOOR) * w
            sev_word = "critical"
        else:
            factor = 1.0 - (1.0 - _INVOLVED_FLOOR) * w
            sev_word = "involved"
        multiplier *= factor
        label = organ.replace("_", "/")
        result["flags"].append(
            f"Serious {label} toxicity conflicts with a {sev_word} "
            f"{label} indication (safety penalty).")

    multiplier = max(_MULTIPLIER_FLOOR, round(multiplier, 3))
    result["multiplier"] = multiplier
    result["penalized"] = multiplier < 0.999
    return result
