"""
Therapeutic context resolver — keeps the platform disease-AGNOSTIC and, crucially,
LOGICALLY COHERENT.

The organ / CNS-relevance is derived from the disease's REAL biology — its Open
Targets **therapeutic areas** (the authoritative disease ontology) — NOT a guess
from the disease name. This is the single source of truth that every dependent
argument (blood-brain-barrier framing, PK framing, route choice) is locked to, so
the platform can never claim e.g. "low brain penetration is a benefit" for a CNS
disease like schizophrenia. A small keyword map is kept only as a last-resort
fallback when Open Targets can't be reached.
"""
import logging
from typing import Dict, List, Optional

from services import http_client

logger = logging.getLogger(__name__)
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# Open Targets therapeutic-area name fragments → canonical organ (order matters;
# CNS is checked first so any nervous-system/psychiatric tag dominates).
_TA_ORGAN = [
    (("nervous system", "psychiatr", "mental health", "cognit", "neurolog", "brain"), "brain"),
    (("eye disease", "ophthalm", "retina"), "eye"),
    (("skin disease", "integument", "cutaneous"), "skin"),
    (("respiratory", "thoracic", "lung", "pulmonary"), "lung"),
    (("cardiovascular", "cardiac", "heart"), "heart"),
    (("digestive", "gastrointestinal", "gastroenter"), "gi"),
    (("musculoskeletal", "connective tissue", "bone disease", "joint"), "joint"),
    (("urinary", "kidney", "renal"), "kidney"),
    (("hematologic", "haematologic", "blood"), "blood"),
    (("reproductive", "breast"), "reproductive"),
    (("endocrine", "metabolic"), "metabolic"),
]

# Last-resort keyword fallback on the raw disease name (only if OT is unavailable).
_KW_ORGAN = [
    ("brain",       ["nervous system", "central nervous", "brain", "neuro", "psychiat", "schizophren",
                     "depress", "bipolar", "psychos", "anxiety", "adhd", "autism", "insomnia",
                     "mental", "cognit", "dementia", "alzheimer", "parkinson", "epilep", "seizure",
                     "migraine", "multiple sclerosis", "huntington"]),
    ("eye",         ["ocular", "retina", "glaucoma", "macular", "cornea", "vision", "ophthalm", "uveitis"]),
    ("skin",        ["skin", "dermat", "psoria", "eczema", "atopic", "acne", "cutaneous", "melanoma"]),
    ("lung",        ["lung", "respiratory", "pulmonary", "asthma", "copd", "bronch", "airway", "cystic fibrosis"]),
    ("heart",       ["cardiac", "cardiovascular", "coronary", "myocard", "arrhythm", "heart failure", "angina"]),
    ("vasculature", ["vascular", "erectile", "hypertension", "hypotension", "thrombo", "arterial", "raynaud"]),
    ("liver",       ["hepat", "liver", "nash", "nafld", "cirrhosis"]),
    ("kidney",      ["renal", "kidney", "nephro", "glomerul"]),
    ("gi",          ["gastro", "intestin", "bowel", "colitis", "crohn", "ulcerative", "ibd", "gastric"]),
    ("joint",       ["arthrit", "rheumat", "synov", "cartilage", "joint", "gout"]),
    ("bone",        ["osteoporos", "bone", "skeletal"]),
    ("blood",       ["leukemia", "leukaemia", "lymphoma", "myeloma", "anemia", "haemato", "hemato"]),
    ("metabolic",   ["diabet", "obesity", "metabolic", "dyslipid", "hyperlip"]),
]

# Organs where CNS exposure is undesirable (peripheral); brain is "required";
# the rest ("systemic"/"metabolic") are "neutral".
_PERIPHERAL = {"eye", "skin", "lung", "heart", "vasculature", "liver", "kidney",
               "gi", "joint", "bone", "blood", "reproductive"}

_ROUTE_HINT = {
    "brain": ["intranasal", "iv", "oral"], "eye": ["ocular", "oral"], "skin": ["topical", "oral"],
    "lung": ["inhalation", "oral"], "heart": ["oral", "iv"], "vasculature": ["oral", "iv"],
    "liver": ["oral", "iv"], "kidney": ["oral", "iv"], "gi": ["oral"], "joint": ["oral", "im"],
    "bone": ["oral", "iv"], "blood": ["iv", "oral"], "reproductive": ["oral"],
    "metabolic": ["oral"], "systemic": ["oral", "iv"],
}

_areas_cache: Dict[str, list] = {}


def _fetch_areas(disease: str) -> List[str]:
    """Open Targets therapeutic areas for a disease (resolve name → EFO → areas)."""
    if not disease:
        return []
    key = disease.strip().lower()
    if key in _areas_cache:
        return _areas_cache[key]
    areas: List[str] = []
    try:
        efo = ""
        try:
            from services.disease_ontology import resolve_disease
            efo = (resolve_disease(disease) or {}).get("disease_id", "")
        except Exception:
            efo = ""
        if efo:
            q = "query($id:String!){ disease(efoId:$id){ therapeuticAreas{ name } } }"
            r = http_client.post(OT_URL, json={"query": q, "variables": {"id": efo}},
                                 timeout=10, headers={"Content-Type": "application/json"})
            if r and r.ok:
                dis = (r.json().get("data", {}) or {}).get("disease") or {}
                areas = [a.get("name", "") for a in (dis.get("therapeuticAreas") or [])]
    except Exception as e:
        logger.debug(f"therapeutic areas fetch failed for {disease}: {e}")
    _areas_cache[key] = areas
    return areas


def _organ_from_areas(areas: List[str]) -> Optional[str]:
    blob = " ".join(areas or []).lower()
    if not blob.strip():
        return None
    for frags, organ in _TA_ORGAN:
        if any(f in blob for f in frags):
            return organ
    return "systemic"          # areas present but unmapped → systemic (not None)


def _organ_from_keywords(disease: str) -> str:
    blob = (disease or "").lower()
    for organ, kws in _KW_ORGAN:
        if any(k in blob for k in kws):
            return organ
    return "systemic"


def therapeutic_context(disease: str = "", therapeutic_areas: Optional[List[str]] = None) -> Dict:
    """Resolve the organ/route context for a disease from its real biology.

    `therapeutic_areas` (e.g. the Open Targets areas already attached to a
    candidate) are used directly when provided — no extra network call. Otherwise
    they are fetched; only if that fails do we fall back to disease-name keywords.
    """
    areas = therapeutic_areas if therapeutic_areas is not None else _fetch_areas(disease)
    organ = _organ_from_areas(areas)
    source = "open_targets"
    if organ in (None, "systemic"):
        kw = _organ_from_keywords(disease)
        if kw != "systemic":            # keyword found something more specific
            organ, source = kw, "keyword_fallback"
        elif organ is None:
            organ, source = "systemic", "keyword_fallback"

    is_cns = organ == "brain"
    if is_cns:
        cns_relevance = "required"      # MUST cross the BBB — penetration is essential
    elif organ in _PERIPHERAL:
        cns_relevance = "liability"     # peripheral — CNS exposure = central side-effects
    else:
        cns_relevance = "neutral"       # systemic/metabolic — don't over-claim either way

    axes = (["cns_mpo", "half_life", "selectivity", "safety"] if is_cns
            else ["half_life", "selectivity", "safety", "oral_bioavailability"])

    return {
        "disease": disease, "organ": organ, "is_cns": is_cns,
        "cns_relevance": cns_relevance, "classification_source": source,
        "therapeutic_areas": areas,
        "route_hint": _ROUTE_HINT.get(organ, ["oral"]),
        "relevant_axes": axes,
        "cns_mpo_label": ("CNS-MPO (BBB)" if is_cns else "CNS exposure"),
        "cns_mpo_higher_is_better": is_cns,
    }
