"""
Reverse Repurposing Engine  —  drug → new indications

The disease→drug engine (services/repurposing_engine.py) answers "which molecules
could treat this disease?". This module answers the inverse: "given a molecule the
company already makes, which NEW indications could it be repurposed into?".

Everything is derived at runtime — there is no hardcoded list of indications per
drug. The candidate diseases come from the molecule's own protein targets:

    drug ──(ChEMBL mechanism)──▶ target genes
         ──(Open Targets target→associatedDiseases)──▶ candidate diseases
         ──(subtract the drug's known indications)──▶ NEW opportunities
         ──(reuse the 6-dimension scorer, inverted)──▶ ranked candidates

An optional therapeutic-area filter (e.g. "dermatology", "ophthalmology") is
applied using Open Targets' own therapeuticAreas annotation on each disease — not
a keyword list maintained in code.
"""

import json
import logging
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional

import requests  # noqa: F401

from services import http_client

logger = logging.getLogger(__name__)

# Card Safety/Efficacy column formatters (display only; no scoring). Fail-soft stubs
# keep the reverse screen working if the miner module is unavailable.
try:
    from services.clinical_evidence import (safety_column as _ct_safety_col,
                                             efficacy_column as _ct_efficacy_col)
except Exception:                                    # pragma: no cover
    def _ct_safety_col(_ae, **_kw):
        return {"status": "no_data", "label": "No trial data", "muted": True}
    def _ct_efficacy_col(_v, **_kw):
        return {"status": "no_data", "label": "No trial data", "muted": True}

OT_URL      = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CT_BASE     = "https://clinicaltrials.gov/api/v2/studies"
NCBI_BASE   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FDA_LABEL_URL = "https://api.fda.gov/drug/label.json"

_fda_label_cache: Dict[str, str] = {}


def _fda_label_text(drug_name: str) -> str:
    """The FDA-approved 'indications and usage' prose for a drug (openFDA label),
    lowercased. This is the AUTHORITATIVE approval source — unlike ChEMBL
    max_phase_for_ind (which is the max TRIAL phase, so an FDA-approved indication whose
    pivotal trial was Phase 3 is wrongly read as 'studied, not approved'). Fail-soft → ''."""
    key = (drug_name or "").strip().lower()
    if not key:
        return ""
    if key in _fda_label_cache:
        return _fda_label_cache[key]
    text = ""
    try:
        for field in ("openfda.generic_name", "openfda.brand_name"):
            r = http_client.get(FDA_LABEL_URL,
                                params={"search": f'{field}:"{key}"', "limit": 1},
                                timeout=10)
            if r and r.ok:
                results = r.json().get("results", [])
                if results:
                    iu = results[0].get("indications_and_usage") or []
                    text = " ".join(iu).lower() if isinstance(iu, list) else str(iu).lower()
                    if text:
                        break
    except Exception as e:
        logger.debug(f"FDA label fetch failed for {drug_name}: {e}")
    _fda_label_cache[key] = text
    return text


def _mark_fda_approved(drug_name: str, known_indications: List[Dict]) -> None:
    """Upgrade an indication's phase to 4 (approved) when the FDA label names it — closes
    the ChEMBL trial-phase gap that mislabels approved biologic indications (mepolizumab→
    EGPA / nasal polyps) as 'studied'. Fail-soft (leaves ChEMBL phase if no label).

    PRECISION (audited 2026-07): disease SYNONYMS (Churg-Strauss = eosinophilic
    granulomatosis with polyangiitis = EGPA) are matched, and a QUALIFIER discriminator
    prevents a false match of the distinct sibling — EGPA (eosinophilic GPA) must NOT
    approve plain GPA/Wegener's just because 'granulomatosis with polyangiitis' is a
    substring of the label's 'EOSINOPHILIC granulomatosis with polyangiitis'."""
    label = _fda_label_text(drug_name)
    if not label:
        return
    # indication name (lower) → extra approved-phrase synonyms to also look for
    _SYN = {"churg-strauss syndrome": "eosinophilic granulomatosis with polyangiitis",
            "churg strauss syndrome": "eosinophilic granulomatosis with polyangiitis"}
    _stop = {"disease", "disorder", "syndrome", "the", "of", "and", "with", "chronic"}
    def _toks(s):
        return {w for w in re.split(r"[^a-z0-9]+", (s or "").lower())
                if len(w) > 3 and w not in _stop}
    ltoks = _toks(label)
    for k in known_indications:
        nm = (k.get("name") or "").lower()
        nt = _toks(nm + " " + _SYN.get(nm, ""))
        if not nt:
            continue
        # 'eosinophilic' qualifier discriminator: if the indication is granulomatosis+
        # polyangiitis WITHOUT 'eosinophilic' (GPA/Wegener) and every 'granulomatosis' in
        # the label is preceded by 'eosinophilic' (i.e. the label only approves EGPA),
        # this is a different disease → do not approve it.
        if ("granulomatosis" in nt and "eosinophilic" not in nt
                and "eosinophilic granulomatosis" in label
                and re.search(r"(?<!eosinophilic )granulomatosis", label) is None):
            continue
        if len(nt & ltoks) / len(nt) >= 0.6:             # indication named in the FDA label
            if float(k.get("max_phase") or 0) < 4:
                k["max_phase"] = 4.0
                k["approval_source"] = "fda_label"
CACHE_FILE  = Path(__file__).parent.parent / "data" / "reverse_cache.json"
CACHE_TTL   = 21600  # 6 hours

# How many trial conditions / literature MeSH diseases to resolve & merge
MAX_TRIAL_CONDITIONS = 20
MAX_LIT_DISEASES      = 15

# How many candidate diseases to fully score (each costs an Open Targets +
# STRING call). Candidates are ranked by association score first, then trimmed.
MAX_SCORED_CANDIDATES = 25
MIN_ASSOCIATION_SCORE = 0.05

# Adaptive "actionable" cutoff (composite score band). The base bar is heuristic
# (pending the RepoDB clean-controls validation); the exceptions are the operational
# logic: a real trial link relaxes it, an orphan disease relaxes it. Exclusion of
# approved indications runs BEFORE ranking, so these can never re-admit an on-label use.
_BASE_ACTIONABLE_CUTOFF = 0.40
_TRIAL_LINK_CUTOFF      = 0.30   # active/historic trial in THIS indication overrides
_ORPHAN_CUTOFF_FACTOR   = 0.85   # rare disease → relax the bar 15%

try:
    from services.disease_ontology import is_orphan_candidate
except Exception:                                    # pragma: no cover
    def is_orphan_candidate(_name: str, _efo: str = "") -> bool:  # fail-soft stub
        return False


# ── Adverse-event / side-effect filtering ────────────────────────────────────
#
# When PubMed indexes an article about "drug X causes/induces symptom Y", the
# MeSH heading for Y carries a QUALIFIER (subheading) marking it as an adverse
# effect — NOT a therapeutic use. Harvesting the bare DescriptorName and treating
# it as a candidate INDICATION is how a drug's own side effects (Headache,
# Nausea, Hypertension, Rash…) were leaking in as diseases to repurpose the drug
# INTO. We inspect each heading's <QualifierName> children to tell the two apart.
#
# These are MeSH subheadings, matched case-insensitively on the QualifierName
# text. This is drug-agnostic reasoning about the SEMANTICS of the co-mention —
# not a hardcoded per-drug indication list — consistent with this module's
# "everything derived at runtime" principle.
ADVERSE_QUALIFIERS = {
    "chemically induced",   # "…, chemically induced" → drug caused the condition
    "adverse effects",
    "toxicity",
    "poisoning",
}
THERAPEUTIC_QUALIFIERS = {
    "drug therapy",
    "therapeutic use",
    "prevention & control",
    "prevention and control",
}

# Curated backstop denylist: MeSH descriptors that are inherently side-effect /
# symptom / lab-finding categories and are essentially NEVER a legitimate
# repurposing indication — regardless of qualifier. This is a GENERIC
# adverse-event / symptom denylist (not drug-specific), so it is consistent with
# the "no hardcoded list of indications per drug" principle. Kept deliberately
# conservative: it lists bare symptom / AE terms, NOT disease entities (we do not
# want to nuke real indications — e.g. "Migraine Disorders" is a disease and is
# intentionally NOT here even though "Headache" is).
AE_SYMPTOM_DENYLIST = {
    "drug-related side effects and adverse reactions",
    "drug hypersensitivity",
    "drug eruptions",
    "nausea",
    "vomiting",
    "headache",
    "dizziness",
    "diarrhea",
    "fatigue",
    "rash",
    "exanthema",          # MeSH descriptor for "rash"
    "pruritus",
    "constipation",
    "edema",
    "oedema",             # British spelling as it appears in some ontologies
    "weight gain",
    "weight loss",
    "somnolence",
    "sleepiness",
    "disorders of excessive somnolence",  # MeSH descriptor for somnolence
    "insomnia",
    "sleep initiation and maintenance disorders",  # MeSH descriptor for insomnia
    # Generic hematologic / organ TOXICITY sign terms. These are the biological
    # CONSEQUENCE of a cytotoxic/kinase-inhibitor mechanism (e.g. CDK4/6i →
    # myelosuppression), never a disease the drug is repurposed to TREAT. Kept
    # deliberately to BARE toxicity signs — real disease ENTITIES that happen to
    # be hematologic (aplastic anemia, myelodysplastic syndrome, acute myeloid
    # leukemia, myelofibrosis, pernicious anemia) are NOT listed and must survive.
    "myelosuppression",
    "bone marrow suppression",
    "neutropenia",
    "febrile neutropenia",
    "leukopenia",
    "thrombocytopenia",
    "pancytopenia",
    "agranulocytosis",
    "anemia",              # bare "anemia" only; specific anemias (aplastic, pernicious…) are diseases
    "hepatotoxicity",
    "nephrotoxicity",
    "cardiotoxicity",
    "qt prolongation",
}

# ── Phenotype-not-disease filtering (Open Targets / ontology) ─────────────────
#
# Open Targets' target→associatedDiseases returns whatever ontology terms are
# associated with a gene, which includes SIGNS / PHENOTYPES / MEASUREMENTS that
# are biological consequences of perturbing the target — not disease entities the
# drug could be repurposed to treat. For palbociclib's targets (CDK4/CDK6/CCND1)
# these are toxicity phenotypes: "myelosuppression" (EFO_0007053), "Decreased
# total neutrophil count" (HP_0001875), "Myelodysplasia" (HP_0002863).
#
# A term is a phenotype-not-disease when EITHER:
#   (a) its ontology id is in the Human Phenotype Ontology (id starts with "HP_")
#       — HPO terms are, by definition, phenotypic abnormalities / signs, or
#   (b) it carries NO therapeutic area, or ALL of its therapeutic areas are
#       non-disease buckets — the generic "phenotype" tag, a "measurement"
#       (e.g. "body weight"), a "biological_process" (e.g. "alcohol drinking"),
#       or "injury, poisoning or other complication" (e.g. "dislocation",
#       bare "injury"). None of these are diseases a drug is repurposed to TREAT.
# A REAL disease that merely ALSO carries one of these tags among several areas
# (e.g. pernicious anemia / myelofibrosis carry "phenotype" but also "hematologic
# disorder"; diabetes carries "phenotype" but also "endocrine system disorder")
# has ≥1 genuine disease area and is correctly KEPT.
_NON_DISEASE_AREAS = {
    "phenotype",
    "measurement",
    "biological_process",
    "injury, poisoning or other complication",
}


def _is_phenotype_not_disease(efo_id: str, therapeutic_areas: List[str]) -> bool:
    """True when an OT/ontology term is a sign/phenotype/measurement/injury, not a disease."""
    if (efo_id or "").upper().startswith("HP_"):
        return True
    areas = [a for a in (therapeutic_areas or []) if (a or "").strip()]
    if not areas:
        return True
    # Drop only when EVERY area is a non-disease bucket (no genuine disease area).
    return all(a.strip().lower() in _NON_DISEASE_AREAS for a in areas)


# Disease identity matching is centralized in services.disease_id (the platform-wide
# reconciliation layer) so ChEMBL MeSH headings, OT MONDO/EFO names, and RepoDB/UMLS
# free-text all match through one implementation. Thin local aliases keep call sites.
from services.disease_id import canonical_key as _disease_tokens, same_disease as _same_disease


def _is_ae_symptom(name: str) -> bool:
    """True if a descriptor/condition name is on the generic AE/symptom denylist."""
    return (name or "").strip().lower() in AE_SYMPTOM_DENYLIST


# ── Cache ───────────────────────────────────────────────────────────────────

# Guards concurrent read-modify-write of the on-disk cache. The reverse/pair
# scorers are now called from worker threads (pathway screen), so an unguarded
# full-file rewrite could interleave and corrupt the JSON or drop sibling entries.
_CACHE_LOCK = threading.Lock()


def _load_cache() -> dict:
    try:
        if CACHE_FILE.exists():
            data = json.loads(CACHE_FILE.read_text(encoding="utf-8"))
            now = time.time()
            return {k: v for k, v in data.items() if now - v.get("_ts", 0) < CACHE_TTL}
    except Exception:
        pass
    return {}


def _save_cache(cache: dict):
    try:
        CACHE_FILE.parent.mkdir(exist_ok=True)
        # Atomic replace so a concurrent reader never sees a half-written file.
        tmp = CACHE_FILE.with_suffix(f".{os.getpid()}.{threading.get_ident()}.tmp")
        tmp.write_text(json.dumps(cache), encoding="utf-8")
        os.replace(tmp, CACHE_FILE)
    except Exception:
        pass


def _cache_put(key: str, value: dict):
    """Append a single entry without clobbering entries written concurrently."""
    with _CACHE_LOCK:
        cache = _load_cache()
        cache[key] = value
        _save_cache(cache)


# ── Open Targets helpers ────────────────────────────────────────────────────

def _ot_graphql(query: str, variables: dict, timeout: int = 12) -> dict:
    try:
        r = http_client.post(OT_URL, json={"query": query, "variables": variables},
                          timeout=timeout, headers={"Content-Type": "application/json"})
        if r and r.ok:
            return r.json().get("data", {})
    except Exception as e:
        logger.debug(f"Open Targets GraphQL error: {e}")
    return {}


@lru_cache(maxsize=4096)
def _gene_to_ensembl(symbol: str) -> str:
    """Resolve a gene symbol → Ensembl target ID via Open Targets search.
    Memoized: the same gene recurs across a drug's targets and across drugs, so caching
    the (stable) symbol→id lookup removes most of the reverse-screen's live-call latency."""
    q = """
    query($q: String!) {
      search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
        hits { id name entity }
      }
    }"""
    hits = _ot_graphql(q, {"q": symbol}).get("search", {}).get("hits", [])
    return hits[0]["id"] if hits else ""


@lru_cache(maxsize=8192)
def _disease_known_drugs(efo_id: str) -> int:
    """Count of drugs + clinical candidates for a disease — a tractability signal.
    Memoized: the same EFO id recurs across candidates and drugs."""
    q = """
    query($id: String!) {
      disease(efoId: $id) { drugAndClinicalCandidates { count } }
    }"""
    data = _ot_graphql(q, {"id": efo_id})
    disease = data.get("disease") or {}
    kd = disease.get("drugAndClinicalCandidates") or {}
    return kd.get("count", 0) or 0


_PHASE_MAP = {"PHASE4": 4, "PHASE3": 3, "PHASE2": 2, "PHASE1": 1, "EARLY_PHASE1": 1}


def _max_phase_from_phases(phases: List[str]) -> int:
    return max((_PHASE_MAP.get((p or "").upper(), 0) for p in phases), default=0)


# Condition strings that are not indications: toxicity/complication reports and
# drug / drug-class labels that trials list as a "condition". These are matched
# generically (substring / suffix), NOT against a per-drug list, so the rule
# holds for any molecule.
_TOXICITY_CONDITION_TOKENS = ("toxicity", "complication", "adverse")


def _is_non_indication_condition(cond: str, drug_name: str = "") -> bool:
    """True when a trial 'condition' string is not a disease to repurpose INTO.

    Drops (a) toxicity / complication / adverse-event conditions, and (b)
    conditions that are actually a drug or drug-CLASS name rather than a disease
    (equal to the queried drug, or a pharmacologic-class label such as
    '... inhibitor' / '... antagonist' / '... agonist' / '... blocker').
    """
    c = (cond or "").strip().lower()
    if not c:
        return True
    if any(tok in c for tok in _TOXICITY_CONDITION_TOKENS):
        return True
    if drug_name and c == drug_name.strip().lower():
        return True
    # Pharmacologic-class / mechanism labels (e.g. "CDK4/6 Inhibitor",
    # "PARP Inhibitor", "Immune Checkpoint Inhibitor") — a drug class, not a disease.
    if re.search(r"\b(inhibitor|antagonist|agonist|blocker|modulator)s?$", c):
        return True
    return False


# P3 trial-OUTCOME classification of a terminated/withdrawn trial's whyStopped
# reason (same taxonomy validated on RepoDB). Efficacy/safety failures are real
# NEGATIVE evidence; business/accrual failures are neutral (not the drug's fault).
_TRIAL_EFF_KW = ("efficac", "futil", "did not meet", "lack of response", "no benefit",
                 "ineffective", "not effective", "endpoint", "lack of effic")
_TRIAL_SAF_KW = ("safety", "adverse", "toxic", "risk/benefit", "risk-benefit",
                 "tolerab", "death", "serious ae")
_TRIAL_BIZ_KW = ("accrual", "enroll", "recruit", "funding", "fund", "sponsor",
                 "business", "strateg", "logist", "supply", "financ", "feasib",
                 "covid", "administrative")


def _classify_why_stopped(why: str) -> str:
    w = (why or "").lower()
    if not w:
        return "unknown"
    if any(k in w for k in _TRIAL_EFF_KW):
        return "efficacy"
    if any(k in w for k in _TRIAL_SAF_KW):
        return "safety"
    if any(k in w for k in _TRIAL_BIZ_KW):
        return "business"
    return "other"


def _trials_for_drug(drug_name: str, page_size: int = 200, probe: dict = None) -> Dict[str, Dict]:
    """Conditions studied in ClinicalTrials.gov for this drug → candidate indications,
    WITH per-condition trial OUTCOMES (P3): completed / with-results / failed-for-
    efficacy / failed-for-safety / ongoing, and a directional outcome_signal
    (-1..+1). A drug that FAILED a trial for efficacy in an indication is real
    negative evidence, not merely "a trial exists".

    Relationship-aware: `query.intr` so the drug is an INTERVENTION, not merely
    mentioned (the v2 param is `query.intr`, not `query.intervention`).

    Stamp-at-fetch: when a mutable `probe` dict is supplied, this helper records
    `probe["ctgov_queried"] = True` only once ClinicalTrials.gov has actually
    ANSWERED (r.ok). That is the honest source of the Verification stage's
    clinicaltrials `checked` flag: an empty answer for a disease with no program
    still reads checked (a query ran), while a network failure leaves it silent.
    """
    out: Dict[str, Dict] = {}
    try:
        r = http_client.get(CT_BASE, params={"query.intr": drug_name,
                                          "pageSize": page_size, "format": "json"}, timeout=15)
        if r and r.ok:
            if probe is not None:
                probe["ctgov_queried"] = True
            for s in r.json().get("studies", []):
                ps = s.get("protocolSection", {})
                conds = ps.get("conditionsModule", {}).get("conditions", []) or []
                phases = ps.get("designModule", {}).get("phases", []) or []
                ph = _max_phase_from_phases(phases)
                st = ps.get("statusModule", {})
                status = (st.get("overallStatus") or "").upper()
                why = st.get("whyStopped", "")
                has_results = bool(s.get("hasResults"))
                stopped = status in ("TERMINATED", "WITHDRAWN", "SUSPENDED")
                fail = _classify_why_stopped(why) if stopped else ""
                for cond in conds:
                    k = cond.strip().lower()
                    if not k or _is_non_indication_condition(cond, drug_name):
                        continue
                    e = out.setdefault(k, {"name": cond.strip(), "count": 0, "max_phase": 0,
                                           "completed": 0, "with_results": 0,
                                           "failed_efficacy": 0, "failed_safety": 0,
                                           "ongoing": 0, "_studies": []})
                    e["count"] += 1
                    e["_studies"].append(s)          # kept for Layer-1 endpoint parsing
                    e["max_phase"] = max(e["max_phase"], ph)
                    if status == "COMPLETED":
                        e["completed"] += 1
                        if has_results:
                            e["with_results"] += 1
                    elif status in ("RECRUITING", "ACTIVE_NOT_RECRUITING",
                                    "ENROLLING_BY_INVITATION", "NOT_YET_RECRUITING"):
                        e["ongoing"] += 1
                    elif fail == "efficacy":
                        e["failed_efficacy"] += 1
                    elif fail == "safety":
                        e["failed_safety"] += 1
            # Layer 1 — Clinical Endpoint Parsing ("did it work?"): classify each
            # condition's trials by their PRIMARY-endpoint outcome (met / biomarker-only /
            # failed / terminated) and let THAT drive the outcome signal — not a raw trial
            # count. Falls back to the honest count-based signal when results aren't parseable.
            try:
                from services.endpoint_parser import aggregate as _ep_agg
            except Exception:
                _ep_agg = None
            try:
                from services.clinical_evidence import parse_adverse_events as _parse_ae
            except Exception:
                _parse_ae = None
            for e in out.values():
                studies = e.pop("_studies", [])
                # Serious-AE (denominatored) summary for the Safety card column — the
                # results-posting study in THIS indication with the largest at-risk
                # denominator. Display only; adds nothing to the composite/rank.
                e["serious_ae"] = None
                if _parse_ae and studies:
                    _best_ae, _best_n = None, -1
                    for _s in studies:
                        _ae = _parse_ae(_s)
                        if _ae is None:
                            continue
                        _m = _ae.get("serious_num_at_risk") or 0
                        if _m > _best_n:
                            _sph = _max_phase_from_phases(
                                (_s.get("protocolSection", {}).get("designModule", {}) or {}).get("phases", []) or [])
                            _best_ae, _best_n = {**_ae, "phase": _sph}, _m
                    e["serious_ae"] = _best_ae
                agg = _ep_agg(studies) if (_ep_agg and studies) else None
                if agg and agg.get("verdict") not in (None, "unknown"):
                    e["outcome_signal"] = agg["outcome_signal"]
                    e["endpoint_verdict"] = agg["verdict"]
                    e["endpoint_note"] = agg.get("note", "")
                    e["endpoint_primary"] = agg.get("primary_endpoint")
                    e["endpoint_p"] = agg.get("p")
                    if agg["verdict"] == "terminated_efficacy":
                        e["failed_efficacy"] = max(e.get("failed_efficacy", 0), 1)
                else:
                    e["outcome_signal"] = _trial_outcome_signal(e)
    except Exception as ex:
        logger.debug(f"trials fetch failed: {ex}")
    return out


def _trial_outcome_signal(e: Dict) -> float:
    """Directional trial-outcome signal in [-1, 1].

    HONESTY GUARD (audited 2026-07): 'has posted results' (ClinicalTrials.gov
    hasResults) means the trial reported SOMETHING — the readout may have been
    NEGATIVE. We do NOT yet parse primary-endpoint significance, so completed /
    with-results / ongoing are credited ONLY as weak PROGRAM signal (a human thought
    the pair worth running), hard-capped so they can never masquerade as a positive
    EFFICACY readout. The only strong signals here are the NEGATIVES — a trial
    stopped for efficacy/safety is real negative evidence."""
    pos = min(0.20, 0.04 * e.get("completed", 0)
                    + 0.05 * e.get("with_results", 0)
                    + 0.03 * e.get("ongoing", 0))
    neg = 0.60 * e.get("failed_efficacy", 0) + 0.40 * e.get("failed_safety", 0)
    return round(max(-1.0, min(1.0, pos - neg)), 3)


def _attach_trials_by_name(c: Dict, trials: Dict) -> None:
    """Attach trial outcomes to a candidate by DISEASE NAME when the EFO-keyed merge
    missed it — a mechanism-generated candidate and its trial often resolve to different
    EFO ids, so a drug that FAILED a trial in this indication would otherwise carry no
    negative. Runs before scoring so the failure actually lowers the composite."""
    if c.get("trial_count") or not trials:
        return
    dn = (c.get("disease") or "").strip().lower()
    if not dn:
        return
    best = trials.get(dn)
    if best is None:
        _norm = lambda s: set(re.sub(r"[^a-z0-9 ]", " ", s.lower()).split()) - {"the", "of", "and"}
        toks = _norm(dn)
        for k, v in trials.items():
            kt = _norm(k)
            if kt and len(kt) >= 2 and (kt <= toks or toks <= kt):   # one name's tokens ⊆ the other
                best = v
                break
    if best:
        c["trial_count"] = best.get("count", 0)
        c["max_trial_phase"] = best.get("max_phase", 0)
        c["trial_outcome_signal"] = best.get("outcome_signal", 0.0)
        c["failed_efficacy"] = best.get("failed_efficacy", 0)
        c["failed_safety"] = best.get("failed_safety", 0)
        c["serious_ae"] = best.get("serious_ae")          # for the Safety card column
        # Layer 1: carry the parsed primary-endpoint verdict onto the candidate
        for _k in ("endpoint_verdict", "endpoint_note", "endpoint_primary", "endpoint_p"):
            if best.get(_k) is not None:
                c[_k] = best.get(_k)
        # Display-only mirror for the Efficacy column (same key the merge path fills)
        if best.get("endpoint_verdict"):
            c["ct_efficacy_verdict"] = best.get("endpoint_verdict")
            c["ct_efficacy_note"] = best.get("endpoint_note", "")
        c["sources"] = set(c.get("sources", set())) | {"clinical_trial"}


def _why_not(e: Dict) -> List[str]:
    """Concise reasons a candidate is down-ranked or non-actionable (the 'why not' view),
    assembled from the gates/penalties already computed for the card."""
    out: List[str] = []
    pkg = e.get("primekg") or {}
    if pkg.get("relation") == "contraindication":
        out.append("Contraindicated: " + (pkg.get("flag")
                   or "PrimeKG labels this an established contraindication for the disease (expert ground truth)."))
    et = e.get("evidence_tier") or {}
    if et.get("tier") == "contradicted":
        out.append("Contradicted: " + (et.get("note") or "a trial failed for efficacy here"))
    appr = e.get("appropriateness") or {}
    if appr.get("appropriate") is False and appr.get("flags"):
        out.append("Wrong direction: " + ", ".join(appr["flags"]))
    saf = e.get("safety") or {}
    if saf.get("penalized") and saf.get("flags"):
        out.append("Safety liability: " + "; ".join(saf["flags"][:2]))
    tf = e.get("trial_failure") or {}
    if tf.get("penalized"):
        out.append("Failed a prior clinical trial in this indication")
    cov = e.get("coverage") or {}
    if cov.get("penalized"):
        out.append(f"Partial target coverage ({cov.get('n_hit', '?')}/{cov.get('n_drivers', '?')} drivers)")
    cc = e.get("clinical_constraints") or {}
    if cc.get("penalized") and cc.get("flags"):
        out.append("Clinical mismatch: " + "; ".join(cc["flags"][:1]))
    cf = e.get("commercial_friction") or {}
    if cf.get("penalized") and cf.get("flags"):
        out.append("Commercial friction: " + cf["flags"][0])
    if (e.get("ctpa") or {}).get("phantom"):
        out.append("No functional cohesion; the target match is off target")
    eb = e.get("evidence_balance") or {}
    if eb.get("verdict") == "insufficient_coverage":
        out.append("Unverified: " + (eb.get("verdict_reason")
                   or "negative sources not yet checked for this pair; hypothesis only."))
    if not e.get("actionable") and not out:
        out.append("Below the evidence threshold: mechanistic signal too weak to act on")
    return out


def _evidence_tier(c: Dict) -> Dict:
    """Honest real-world-evidence strength for a candidate indication, from its trial
    outcomes — so a Phase-III benefit, an unproven early signal, and a failed trial are
    never shown as the same thing. tier ∈ {contradicted, trial-supported, promising,
    literature, preclinical}."""
    fe = c.get("failed_efficacy", 0) or 0
    fs = c.get("failed_safety", 0) or 0
    ph = int(c.get("max_trial_phase", 0) or 0)
    tc = c.get("trial_count", 0) or 0
    lit = c.get("lit_count", 0) or 0
    osig = float(c.get("trial_outcome_signal", 0.0) or 0.0)
    # Layer 1 — parsed PRIMARY-endpoint verdict (highest precision, grounded in the posted
    # p-value / status). When present it decides the tier; else fall through to the
    # honest count-based logic below.
    ev = c.get("endpoint_verdict")
    epn = c.get("endpoint_note", "")
    if ev == "met_primary":
        return {"tier": "trial-supported", "label": "Met primary endpoint",
                "note": epn or "Primary endpoint met in a posted trial."}
    if ev == "terminated_efficacy":
        return {"tier": "contradicted", "label": "Terminated for efficacy",
                "note": epn or "A trial was stopped for lack of efficacy here."}
    if ev == "failed_primary":
        return {"tier": "failed-endpoint", "label": "Failed primary endpoint",
                "note": epn or "The primary clinical endpoint was missed here."}
    if ev == "biomarker_only":
        return {"tier": "biomarker-only", "label": "Biomarker-positive · clinically unproven",
                "note": epn or "A biomarker moved but the primary clinical endpoint was missed."}
    if fe > 0 or osig <= -0.30:
        return {"tier": "contradicted", "label": "Contradicted here",
                "note": "a trial failed for efficacy in this indication"}
    if fs > 0:
        return {"tier": "contradicted", "label": "Failed on safety",
                "note": "a trial stopped for safety in this indication"}
    # A completed/with-results trial we could NOT parse is prior clinical art (a human
    # tested it), NOT proof of benefit — label it exactly that.
    if ph >= 2 and tc > 0:
        return {"tier": "tested-unverified", "label": "Clinically tested · unproven",
                "note": (f"Phase {ph} trial(s) run in this indication, but the primary-endpoint "
                         f"outcome is not verified in-platform — prior clinical art, not proof of benefit.")}
    if tc > 0 or ph >= 1:
        return {"tier": "promising", "label": "Early clinical · unproven",
                "note": f"{tc} trial(s) up to Phase {ph}; no verified positive readout."}
    if lit >= 3:
        return {"tier": "literature", "label": "Literature signal",
                "note": f"{lit} co-mentions, no clinical trial"}
    return {"tier": "preclinical", "label": "Mechanistic only",
            "note": "genetic / pathway association only"}


def _literature_for_drug(drug_name: str, retmax: int = 120) -> Dict[str, int]:
    """Disease MeSH headings co-occurring with the drug in PubMed → candidate indications."""
    counts: Dict[str, int] = {}
    try:
        sr = http_client.get(f"{NCBI_BASE}/esearch.fcgi",
                          params={"db": "pubmed", "term": f'"{drug_name}"[tiab]',
                                  "retmax": retmax, "retmode": "json"}, timeout=10)
        ids = (sr.json() if sr and sr.ok else {}).get("esearchresult", {}).get("idlist", [])
        if not ids:
            return counts
        fr = http_client.get(f"{NCBI_BASE}/efetch.fcgi",
                          params={"db": "pubmed", "id": ",".join(ids), "retmode": "xml"},
                          timeout=25)
        if fr and fr.ok:
            import xml.etree.ElementTree as ET
            root = ET.fromstring(fr.content)
            for mh in root.iter("MeshHeading"):
                d = mh.find("DescriptorName")
                if d is None or not (d.text or "").strip():
                    continue
                name = d.text.strip()

                # Generic AE/symptom backstop: bare symptom / adverse-reaction
                # descriptors are never repurposing indications, drop regardless
                # of qualifier.
                if _is_ae_symptom(name):
                    continue

                # Relationship-aware, NOT keyword co-occurrence: literature
                # contributes a candidate ONLY via an explicit drug→disease
                # THERAPEUTIC relation. A disease MeSH seeds a candidate indication
                # solely when this article frames it as something the drug is used
                # to TREAT — i.e. the heading carries a therapeutic qualifier (drug
                # therapy / therapeutic use / prevention & control). A bare
                # co-mention with NO qualifier is not evidence of an indication and
                # is dropped — this alone removes the Humans/Female/Aged demographic
                # and Piperazines/Protein-Kinase-Inhibitors chemical-class noise,
                # since none of those ever carry a therapeutic qualifier. Adverse-
                # qualified headings (chemically induced / toxicity / poisoning /
                # adverse effects) are likewise excluded.
                quals = {(q.text or "").strip().lower()
                         for q in mh.findall("QualifierName") if (q.text or "").strip()}
                has_therapeutic = bool(quals & THERAPEUTIC_QUALIFIERS)
                has_adverse = bool(quals & ADVERSE_QUALIFIERS)
                if not has_therapeutic or has_adverse:
                    continue

                counts[name] = counts.get(name, 0) + 1
    except Exception as e:
        logger.debug(f"literature fetch failed: {e}")
    return counts


def _resolve_disease_meta(name: str) -> Dict:
    """Resolve a free-text condition / MeSH term → {efo, name, areas} via Open Targets search."""
    q = """
    query($q: String!) {
      search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 1}) {
        hits { id name object { __typename ... on Disease { therapeuticAreas { name } } } }
      }
    }"""
    hits = _ot_graphql(q, {"q": name}).get("search", {}).get("hits", [])
    if not hits:
        return {}
    h = hits[0]
    obj = h.get("object") or {}
    if obj.get("__typename") != "Disease":
        return {}
    return {
        "efo":   h.get("id", ""),
        "name":  h.get("name", name),
        "areas": [t.get("name", "") for t in (obj.get("therapeuticAreas") or [])],
    }


@lru_cache(maxsize=4096)
def _diseases_for_target(ensembl_id: str, size: int = 250) -> List[Dict]:
    """Open Targets diseases associated with a target, with therapeutic areas.
    Memoized (callers read, never mutate the returned rows)."""
    q = """
    query($id: String!, $size: Int!) {
      target(ensemblId: $id) {
        associatedDiseases(page: {index: 0, size: $size}) {
          rows {
            score
            disease { id name therapeuticAreas { id name } }
          }
        }
      }
    }"""
    data = _ot_graphql(q, {"id": ensembl_id, "size": size})
    target = data.get("target") or {}
    rows = (target.get("associatedDiseases") or {}).get("rows", [])
    out = []
    for r in rows:
        d = r.get("disease") or {}
        if not d.get("id"):
            continue
        out.append({
            "efo_id":            d["id"],
            "name":              d.get("name", ""),
            "therapeutic_areas": [ta.get("name", "") for ta in (d.get("therapeuticAreas") or [])],
            "score":             round(r.get("score", 0.0), 4),
        })
    return out


# ── ChEMBL helpers ──────────────────────────────────────────────────────────

def _molecule_detail(chembl_id: str) -> Dict:
    try:
        r = http_client.get(f"{CHEMBL_BASE}/molecule/{chembl_id}.json", timeout=10)
        if r and r.ok:
            return r.json()
    except Exception:
        pass
    return {}


def _to_parent(chembl_id: str) -> str:
    """Map a salt/child molecule to its parent — salts often carry no mechanism data."""
    m = _molecule_detail(chembl_id)
    parent = (m.get("molecule_hierarchy") or {}).get("parent_chembl_id")
    return parent or chembl_id


def _molecule_forms(parent_id: str) -> List[str]:
    """All molecule forms (parent + salts) for a parent — ChEMBL records mechanisms
    inconsistently across forms, so we must look at every form to find the targets."""
    forms = [parent_id]
    try:
        r = http_client.get(f"{CHEMBL_BASE}/molecule_form.json",
                         params={"parent_chembl_id": parent_id, "limit": 50, "format": "json"},
                         timeout=10)
        if r and r.ok:
            for f in r.json().get("molecule_forms", []):
                mid = f.get("molecule_chembl_id")
                if mid and mid not in forms:
                    forms.append(mid)
    except Exception:
        pass
    return forms


def _local_indications(chembl_ids: List[str]) -> List[Dict]:
    """Known indications from the local chembl_33 drug_indication table (salt→parent),
    as [{name, efo_id, max_phase}] — the best phase per indication label. Primary source:
    it's more complete and reliable than the flaky live drug_indication API (e.g. it has
    nintedanib's IPF at phase 4, which the API snapshot misses). [] if the DB is absent."""
    try:
        from services.repurposing_engine import _get_chembl_pool
    except Exception:
        return []
    pool = _get_chembl_pool()
    ids = [c for c in chembl_ids if c]
    if pool is None or not ids:
        return []
    conn = None
    try:
        conn = pool.getconn()
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT DISTINCT COALESCE(mh.parent_molregno, md.molregno), md.molregno
                FROM molecule_dictionary md
                LEFT JOIN molecule_hierarchy mh ON mh.molregno = md.molregno
                WHERE md.chembl_id = ANY(%s)
                """, (ids,))
            molregnos = set()
            for parent, mol in cur.fetchall():
                molregnos.update(m for m in (parent, mol) if m is not None)
            if not molregnos:
                return []
            cur.execute(
                """
                SELECT di.mesh_heading, di.efo_id, MAX(di.max_phase_for_ind)
                FROM drug_indication di
                WHERE di.molregno = ANY(%s) AND di.mesh_heading IS NOT NULL
                GROUP BY di.mesh_heading, di.efo_id
                """, (list(molregnos),))
            by_label: Dict[str, Dict] = {}
            for label, efo, ph in cur.fetchall():
                if not label:
                    continue
                key = label.lower()
                try:
                    ph = float(ph or 0)
                except (TypeError, ValueError):
                    ph = 0.0
                efo_id = (efo or "").replace(":", "_")
                if key in by_label:
                    by_label[key]["max_phase"] = max(by_label[key]["max_phase"], ph)
                else:
                    by_label[key] = {"name": label, "efo_id": efo_id, "max_phase": ph}
            return list(by_label.values())
    except Exception as e:
        logger.debug(f"local indications failed: {e}")
        return []
    finally:
        if conn is not None:
            pool.putconn(conn)


def _local_chembl_id_for_name(name: str) -> str:
    """Resolve a drug NAME to a ChEMBL ID from the local chembl_33 DB (pref_name, then
    synonym), preferring the highest development phase. Primary path so a flaky live
    ChEMBL search API can't block resolution. '' if not found or no local DB."""
    try:
        from services.repurposing_engine import _get_chembl_pool
    except Exception:
        return ""
    pool = _get_chembl_pool()
    q = (name or "").strip().lower()
    if pool is None or not q:
        return ""
    conn = None
    try:
        conn = pool.getconn()
        with conn.cursor() as cur:
            cur.execute("SELECT chembl_id FROM molecule_dictionary WHERE LOWER(pref_name) = %s "
                        "ORDER BY max_phase DESC NULLS LAST LIMIT 1", (q,))
            row = cur.fetchone()
            if row and row[0]:
                return row[0]
            cur.execute("SELECT md.chembl_id FROM molecule_synonyms ms "
                        "JOIN molecule_dictionary md ON md.molregno = ms.molregno "
                        "WHERE LOWER(ms.synonyms) = %s ORDER BY md.max_phase DESC NULLS LAST LIMIT 1", (q,))
            row = cur.fetchone()
            return row[0] if row and row[0] else ""
    except Exception as e:
        logger.debug(f"local chembl name resolve failed: {e}")
        return ""
    finally:
        if conn is not None:
            pool.putconn(conn)


def _local_molecule_bundle(chembl_id: str) -> Optional[Dict]:
    """Everything resolve_drug needs from the local chembl_33 DB for a ChEMBL id — parent
    chembl_id, all salt/parent forms, display name, max_phase, SMILES — in one lookup.
    None if the molecule isn't in the local DB (caller falls back to the live ChEMBL API)."""
    try:
        from services.repurposing_engine import _get_chembl_pool
    except Exception:
        return None
    pool = _get_chembl_pool()
    cid = (chembl_id or "").strip().upper()
    if pool is None or not cid:
        return None
    conn = None
    try:
        conn = pool.getconn()
        with conn.cursor() as cur:
            cur.execute("SELECT md.molregno, COALESCE(mh.parent_molregno, md.molregno) "
                        "FROM molecule_dictionary md "
                        "LEFT JOIN molecule_hierarchy mh ON mh.molregno = md.molregno "
                        "WHERE md.chembl_id = %s", (cid,))
            row = cur.fetchone()
            if not row:
                return None
            molregno, parent_molregno = row[0], row[1]
            cur.execute("SELECT chembl_id, pref_name, max_phase FROM molecule_dictionary "
                        "WHERE molregno = %s", (parent_molregno,))
            prow = cur.fetchone()
            parent_chembl = (prow[0] if prow and prow[0] else cid)
            name = (prow[1] if prow and prow[1] else "")
            try:
                max_phase = float(prow[2]) if prow and prow[2] is not None else 0.0
            except (TypeError, ValueError):
                max_phase = 0.0
            cur.execute("SELECT md.chembl_id FROM molecule_dictionary md "
                        "LEFT JOIN molecule_hierarchy mh ON mh.molregno = md.molregno "
                        "WHERE md.molregno = %s OR mh.parent_molregno = %s",
                        (parent_molregno, parent_molregno))
            forms = [r[0] for r in cur.fetchall() if r[0]]
            for x in (cid, parent_chembl):
                if x not in forms:
                    forms.insert(0, x)
            smiles = ""
            for mr in (parent_molregno, molregno):
                cur.execute("SELECT canonical_smiles FROM compound_structures WHERE molregno = %s", (mr,))
                srow = cur.fetchone()
                if srow and srow[0]:
                    smiles = srow[0]
                    break
            return {"parent_chembl_id": parent_chembl, "forms": forms,
                    "name": name, "max_phase": max_phase, "smiles": smiles}
    except Exception as e:
        logger.debug(f"local molecule bundle failed: {e}")
        return None
    finally:
        if conn is not None:
            pool.putconn(conn)


def resolve_drug(drug: str) -> Dict:
    """
    Resolve a drug name or ChEMBL ID to {chembl_id, name, max_phase, targets, known_indications}.
    Fully LOCAL-FIRST over the chembl_33 DB — name→ID, salt→parent, forms, name/phase/SMILES,
    targets, and indications all read from :5433. The live ChEMBL web API (which intermittently
    500s) is only a fallback, so EBI's uptime can no longer block a resolution.
    """
    chembl_id = ""

    if re.fullmatch(r"CHEMBL\d+", drug.strip(), re.I):
        chembl_id = drug.strip().upper()
    else:
        # Local chembl_33 first — the live ChEMBL search API is intermittently 500-ing,
        # which was blocking resolution of drugs the local DB has (e.g. nintedanib).
        chembl_id = _local_chembl_id_for_name(drug)
    if not chembl_id and not re.fullmatch(r"CHEMBL\d+", drug.strip(), re.I):
        try:
            r = http_client.get(f"{CHEMBL_BASE}/molecule/search.json",
                             params={"q": drug, "limit": 25}, timeout=15)
            if r and r.ok:
                mols = r.json().get("molecules", [])
                # The first hit isn't always the canonical drug (e.g. a polyglutamate
                # or salt). Prefer an exact name match with the highest development phase.
                ql = drug.strip().lower()
                def _phase_num(m):
                    try:
                        return float(m.get("max_phase") or 0)
                    except (TypeError, ValueError):
                        return 0.0
                def _pref(m):
                    pref = (m.get("pref_name") or "").lower()
                    return (1 if pref == ql else 0, _phase_num(m))
                if mols:
                    chembl_id = max(mols, key=_pref).get("molecule_chembl_id", "")
        except Exception:
            pass

    if not chembl_id:
        return {"chembl_id": "", "name": drug, "max_phase": 0,
                "targets": [], "known_indications": []}

    # Parent + salt forms + name/phase/SMILES: local chembl_33 first (one query), live
    # ChEMBL API only if the molecule isn't in the local DB. ChEMBL records mechanisms /
    # indications inconsistently across forms, so we keep every form.
    bundle = _local_molecule_bundle(chembl_id)
    if bundle:
        parent_id = bundle["parent_chembl_id"] or chembl_id
        ids = list(bundle["forms"]) or [chembl_id]
        name = bundle["name"] or drug
        max_phase = bundle["max_phase"] or 0
        smiles = bundle["smiles"] or ""
    else:
        parent_id = _to_parent(chembl_id)
        ids = _molecule_forms(parent_id)
        name, max_phase, smiles = drug, 0, ""
        for cid in (parent_id, chembl_id):
            d = _molecule_detail(cid)
            smiles = smiles or (d.get("molecule_structures") or {}).get("canonical_smiles", "")
            if d.get("pref_name"):
                name = d["pref_name"]
                max_phase = d.get("max_phase") or 0
                break
    if chembl_id not in ids:
        ids.insert(0, chembl_id)

    targets: List[str] = []
    try:
        from services.repurposing_engine import _chembl_targets_for_molecules
        gene_map = _chembl_targets_for_molecules(ids)
        for cid in ids:
            for g in gene_map.get(cid, []):
                if g not in targets:
                    targets.append(g)
    except Exception as e:
        logger.debug(f"target fetch failed: {e}")

    # Carry max_phase per indication — novelty must distinguish an APPROVED use
    # (phase 4, blocks a same-use claim) from one merely STUDIED at a low phase
    # (e.g. a drug-interaction/PK co-mention), which does not. Keep the best phase
    # seen per indication. Consumed by services.regulatory_verdict.
    # Primary: local chembl_33 (complete + reliable). Fallback: live API only if local
    # gives nothing — so a flaky API call no longer blanks known_indications.
    known_indications: List[Dict] = _local_indications(ids)
    if not known_indications:
        by_label: Dict[str, Dict] = {}
        for cid in ids:
            try:
                r = http_client.get(f"{CHEMBL_BASE}/drug_indication.json",
                                 params={"molecule_chembl_id": cid, "limit": 100,
                                         "format": "json"}, timeout=10)
                if r and r.ok:
                    for ind in r.json().get("drug_indications", []):
                        label = ind.get("mesh_heading") or ind.get("efo_term") or ""
                        if not label:
                            continue
                        efo = (ind.get("efo_id") or "").replace(":", "_")
                        try:
                            ph = float(ind.get("max_phase_for_ind") or 0)
                        except (TypeError, ValueError):
                            ph = 0.0
                        key = label.lower()
                        if key in by_label:
                            by_label[key]["max_phase"] = max(by_label[key]["max_phase"], ph)
                        else:
                            by_label[key] = {"name": label, "efo_id": efo, "max_phase": ph}
            except Exception:
                pass
        known_indications = list(by_label.values())

    # Correct approval status against the authoritative FDA label (ChEMBL per-indication
    # phase is trial-phase, not approval — so it mislabels approved biologic indications).
    try:
        _mark_fda_approved(name, known_indications)
    except Exception as e:
        logger.debug(f"FDA approval mark skipped: {e}")

    return {
        "chembl_id":         parent_id,
        "name":              name,
        "max_phase":         max_phase,
        "smiles":            smiles,
        "targets":           targets,
        "known_indications": known_indications,
    }


def indication_phase_for(disease: str, known_indications) -> Optional[int]:
    """Development phase of the drug's known indication that BEST matches `disease`
    (token overlap), or None if none matches. Used to PHASE-SCALE prior clinical art in
    scoring so a merely-STUDIED indication (imatinib→pulmonary hypertension, Phase 2/3)
    is not credited like an APPROVED one (imatinib→chronic myeloid leukemia, Phase 4).

    Best-match (not max over all loose matches) so an approved form is not credited to a
    merely-studied sibling: 'acute myeloid leukemia' matches its own ChEMBL entry, not the
    Phase-4 'chronic myeloid leukemia'. 'chronic'/'acute'/'primary' are kept as tokens
    precisely because they discriminate disease subtypes."""
    _stop = {"disease", "disorder", "syndrome", "with", "associated", "the", "of", "and"}
    def _tok(s):
        # keep numeric tokens (len 1) so subtype discriminators survive: 'Type 2' vs
        # 'Type 1' must NOT collapse (that gave T2D the Phase-3 Type-1 entry, breaking
        # metformin→T2D's own-therapy suppression → false 'Contraindicated').
        return {w for w in re.split(r"[^a-z0-9]+", (s or "").lower())
                if (len(w) > 3 or w.isdigit()) and w not in _stop}
    dt = _tok(disease)
    if not dt:
        return None
    best_ph, best_ov = None, -1.0
    for k in known_indications or []:
        nm = k.get("name", "") if isinstance(k, dict) else str(k)
        kt = _tok(nm)
        if not kt:
            continue
        ov = len(dt & kt) / len(dt)          # fraction of disease tokens present in the indication
        _ph = int(float(k.get("max_phase") or 0)) if isinstance(k, dict) else 0
        # higher overlap wins; on a tie prefer the higher (more-approved) phase, so an
        # approved form is never credited to a merely-studied same-name sibling.
        if ov >= 0.5 and (ov > best_ov or (ov == best_ov and _ph > (best_ph or -1))):
            best_ov, best_ph = ov, _ph
    return best_ph


def narrow_broad_disease(disease: str, drug_genes, _cache: Dict = {}) -> Optional[str]:
    """If `disease` is a BROAD umbrella (vasculitis), return the single drug-RELEVANT
    subtype — the one whose disease genes include a drug target (EGPA for an IL-5 drug) —
    so a candidate names the EXACT disease, not the umbrella (mepolizumab does not treat
    vasculitis broadly, only eosinophilic/EGPA). None if not broad or no target-linked
    subtype. Memoized; fail-soft."""
    dg = {g.upper() for g in (drug_genes or [])}
    if not dg or not disease:
        return None
    ck = (disease.lower(), tuple(sorted(dg)))
    if ck in _cache:
        return _cache[ck]
    res = None
    try:
        from services.disease_normalize import expand
        from services.disease_ontology import resolve_disease as _rd
        e = expand(disease)
        if e.get("is_broad"):
            best, best_n = None, 0
            for sub in (e.get("indications") or [])[:6]:
                lbl = sub.get("label", "")
                if not lbl:
                    continue
                sg = {(t.get("gene_symbol") or "").upper()
                      for t in (_rd(lbl) or {}).get("targets", [])[:40]}
                n = len(dg & sg)
                if n > best_n:
                    best_n, best = n, lbl
            res = best
    except Exception as ex:
        logger.debug(f"narrow_broad_disease({disease}) skipped: {ex}")
    _cache[ck] = res
    return res


def _reverse_provenance() -> Dict:
    """Data lineage for the drug→indications ranking: the sources that produced it, each
    dated + Data-Age & Integrity scored. Fail-soft → {}."""
    try:
        from services.provenance import lineage
        return lineage(["drugs_fda_label", "clinicaltrials", "open_targets",
                        "chembl33", "drkg", "mondo"])
    except Exception:
        return {}


# ── Canonical repurposing score (single source of truth for a drug–disease pair) ──

def canonical_pair_score(chembl_id: str, disease: str, drug_genes: Optional[List[str]] = None,
                         max_phase: Optional[int] = None, indications: str = "",
                         drug_name: str = "", mechanistic_prior: Optional[float] = None) -> Dict:
    """THE repurposing score for a (drug, disease) pair — used identically by the
    Repurpose card, the analysis header, and the Evidence Dossier so the number is
    the same everywhere. Composite of the 6 screens; cached per pair (deterministic).

    mechanistic_prior (0..1): optional externally-established drug→disease mechanistic
    linkage from the surfacing engine (pathway modulation, novel-target hit). Folded
    into the pathway dimension so a candidate is not scored as if it had zero mechanism
    when the very reason it surfaced IS a mechanism. Changes the cache key so a
    prior-boosted score never overwrites the plain direct-scoring value for the pair."""
    _pk = "" if mechanistic_prior is None else f"|m{round(float(mechanistic_prior), 2)}"
    # v3: FIX 1 argument-parity (adds therapeutic_areas + studied_for_disease to the canonical
    # scoring call so the dossier matches the forward list) — bump invalidates pre-fix cached
    # pair scores. v2 was phase-scaled prior-art + mechanism-only renorm (2026-07 audit).
    cache_key = f"pair:v3:{(chembl_id or drug_name).lower()}|{disease.lower().strip()}{_pk}"
    cache = _load_cache()
    if cache_key in cache:
        return cache[cache_key]

    try:
        from services.repurposing_engine import score_compound_for_disease
        from services.disease_ontology import resolve_disease, get_ppi_network
    except Exception as e:
        return {"score": 0.0, "scores": {}, "composite_score": 0.0, "error": str(e)}

    # Drug context — fetch once, the same way the reverse engine does, so inputs match
    if drug_genes is None or max_phase is None or not indications:
        info = resolve_drug(chembl_id or drug_name)
        if drug_genes is None:
            drug_genes = info.get("targets", [])
        if max_phase is None:
            max_phase = info.get("max_phase", 0)
        if not indications:
            indications = "; ".join(k["name"] for k in info.get("known_indications", []))
        if not drug_name:
            drug_name = info.get("name", "")

    dinfo = resolve_disease(disease) or {}
    disease_genes    = [t["gene_symbol"] for t in dinfo.get("targets", [])[:40]]
    disease_pathways = dinfo.get("pathways", [])
    ppi = get_ppi_network(disease_genes[:15]) if disease_genes else {}
    compound = {"chembl_id": chembl_id, "name": drug_name, "max_phase": max_phase or 0,
                "indications": indications, "targets": ";".join(drug_genes)}
    # Mechanism action types → direction-aware pathway score (kept consistent with
    # the discovery screen so the same pair scores identically everywhere).
    drug_actions = None
    if chembl_id:
        try:
            from services.repurposing_engine import _actions_for_molecules
            drug_actions = _actions_for_molecules([chembl_id]).get(chembl_id)
        except Exception:
            drug_actions = None
    # Genetic-association weights (drive both the weighted target overlap and the
    # target-coverage gate) from the resolved disease targets.
    disease_weights = {t["gene_symbol"].upper():
                       (t.get("quality_score") or t.get("genetic_score") or t.get("score", 0.0))
                       for t in dinfo.get("targets", []) if t.get("gene_symbol")}
    # per-indication phase (drug's development phase in THIS disease) → phase-scale prior art
    _iphase = None
    try:
        _pinfo = resolve_drug(chembl_id or drug_name)
        _iphase = indication_phase_for(disease, _pinfo.get("known_indications", []))
    except Exception:
        pass
    # FIX 1 — argument PARITY with the forward screen (run_repurposing_screen) so the SAME
    # (drug, disease) pair yields the SAME composite on the leads list and in the dossier.
    # therapeutic_areas revives the LoF/congenital appropriateness gate here (was []); it is
    # resolved from the SAME OT disease-area helper the forward path uses.
    _disease_areas: List[str] = []
    try:
        from services.disease_ontology import therapeutic_areas as _ot_areas
        _disease_areas = _ot_areas(dinfo.get("disease_id", "")) or []
    except Exception:
        _disease_areas = []
    # studied_for_disease: any development footprint (phase>=1) for THIS disease — the SAME
    # confounding-by-indication signal the forward screen passes (approved_chembls_for_disease,
    # min_phase=1). Batched membership query; fail-soft to False.
    _studied = False
    try:
        from services.repurposing_scorer import approved_chembls_for_disease
        if chembl_id:
            _studied = chembl_id in approved_chembls_for_disease([chembl_id], [], disease, min_phase=1)
    except Exception:
        _studied = False
    sr = score_compound_for_disease(compound, disease, disease_genes, disease_pathways, ppi,
                                    drug_genes, drug_actions=drug_actions,
                                    disease_gene_weights=disease_weights or None,
                                    mechanistic_prior=mechanistic_prior,
                                    indication_phase=_iphase,
                                    therapeutic_areas=_disease_areas,
                                    studied_for_disease=_studied)
    out = {"score": sr["composite_score"], "composite_score": sr["composite_score"],
           "scores": sr["scores"], "safety": sr.get("safety", {}),
           "coverage": sr.get("coverage", {}),
           "clinical_constraints": sr.get("clinical_constraints", {}),
           "trial_failure": sr.get("trial_failure", {}),
           "directional_evidence": sr.get("directional_evidence", {}),
           "proliferation": sr.get("proliferation", {}),
           "signature_reversal": sr.get("signature_reversal", {}),
           "direction": sr.get("direction", {}),
           "calibration": sr.get("calibration", {}),
           "primekg": sr.get("primekg", {}),   # ground-truth relation + P(treats) cross-check
           "drug_genes": drug_genes[:20], "_ts": time.time()}

    # Validated mechanistic-plausibility probability (DWPC/metapath model; held-out
    # AUC 0.98 vs random). Additive + fail-soft: present only where the pair maps
    # into the Hetionet subgraph, and clearly scoped as plausibility/triage — NOT a
    # probability of clinical success (which is not learnable from the graph).
    try:
        from services.repurposing_predictor import plausibility
        pl = plausibility(drug_name, disease)
        if pl:
            out["plausibility"] = pl
    except Exception as e:
        logger.debug(f"plausibility skipped: {e}")

    # Modern-KG plausibility (DRKG DistMult) — covers ~5.1k diseases / newer drugs the
    # 2016 Hetionet model misses. Kept as its own axis; also promoted to the primary
    # `plausibility` field when Hetionet has no coverage for this pair, so newer drugs
    # (Palbociclib etc.) still get a KG plausibility instead of a blank.
    try:
        from services.drkg_predictor import treat_probability
        kgp = treat_probability(drug=drug_name, disease=disease, chembl_id=chembl_id)
        if kgp is not None:
            out["kg_plausibility"] = kgp
            if not out.get("plausibility"):
                out["plausibility"] = kgp
    except Exception as e:
        logger.debug(f"DRKG plausibility skipped: {e}")
    _cache_put(cache_key, out)
    return out


def _plausibility_for(drug_name: str, disease: str):
    """Validated DWPC plausibility P(treats) for a pair, or None off Hetionet
    coverage. Fail-soft wrapper so the reverse screen never breaks on it."""
    try:
        from services.repurposing_predictor import plausibility
        return plausibility(drug_name, disease)
    except Exception:
        return None


def _disease_value_for(disease_name: str, efo_id: str = ""):
    """Repurposing Value Score for the candidate INDICATION — 'would a pharma company
    pursue this disease at all' (burden × unmet-need × market-fit). Fail-soft → None."""
    try:
        from services.disease_value import value_for
        return value_for(disease_name, efo_id)
    except Exception:
        return None


def _lead_viability_for(chembl_id: str, name: str, disease: str, genes, smiles: str = "",
                        plaus_p=None):
    """Lead-viability funnel (potency gate here; PBPK window is deferred to the
    dossier for speed). Fail-soft — never breaks the screen."""
    try:
        from services.lead_viability import assess
        return assess(name, chembl_id, disease, list(genes or []), smiles=smiles,
                      plausibility_p=plaus_p, run_window=False)
    except Exception:
        return None


# ── Candidate generation + scoring ──────────────────────────────────────────

# Ranking weights for the REVERSE direction. The drug's GLOBAL clinical/regulatory
# status is deliberately excluded — it's constant for a fixed drug and doesn't
# discriminate between indications. What we rank on is (a) mechanistic evidence —
# shared targets, pathways, network proximity, OT genetic association — and (b)
# real-world repurposing signal SPECIFIC to each candidate indication: existing
# clinical trials of this drug in that disease and literature co-mention. Unlike
# global phase, that clinical signal varies per indication, so it belongs here.
EVIDENCE_WEIGHTS = {
    "target":      0.30,   # Jaccard(drug targets, disease genes)
    "association": 0.20,   # Open Targets target→disease association score
    "pathway":     0.12,
    "ppi":         0.08,
    "clinical":    0.30,   # trials of THIS drug in THIS disease + literature co-mention
}


def _matches_area(therapeutic_areas: List[str], area_filter: str) -> bool:
    """Dynamic area match using Open Targets' own therapeuticAreas annotation."""
    if not area_filter:
        return True
    af = area_filter.lower().strip()
    # Accept a few common synonyms without hardcoding disease lists. Oncology is the
    # important one: Open Targets labels cancers "cancer or benign tumor", NOT "oncology",
    # so the bare filter matched nothing (0 candidates → "no oncology indication").
    synonyms = {
        "dermatology": ["skin", "dermat", "integument"],
        "ophthalmology": ["eye", "ophthalm", "visual", "ocular"],
        "oncology": ["cancer or benign tumor", "cancer", "neoplasm", "tumor", "tumour",
                     "carcinoma", "leukemia", "leukaemia", "lymphoma", "sarcoma",
                     "melanoma", "myeloma", "oncolog"],
    }
    needles = [af] + synonyms.get(af, [])
    blob = " ".join(therapeutic_areas).lower()
    return any(n in blob for n in needles)


def _is_oncology_area(therapeutic_areas: List[str]) -> bool:
    blob = " ".join(therapeutic_areas).lower()
    return "cancer or benign tumor" in blob or "neoplasm" in blob


def _area_is_oncology(area_filter: str) -> bool:
    af = (area_filter or "").lower()
    return "oncolog" in af or "cancer" in af or "neoplasm" in af


# Disease-name oncology check — used to tag DRKG candidates and to classify a drug's
# own indications (a MeSH heading / disease name, not an OT therapeutic-area label).
_ONCOLOGY_TERMS = ("neoplasm", "carcinoma", "cancer", "leukemia", "lymphoma", "melanoma",
                   "sarcoma", "tumor", "tumour", "glioma", "glioblastoma", "myeloma",
                   "blastoma", "mesothelioma", "adenocarcinoma", "malignan")


def _is_oncology_name(name: str) -> bool:
    n = (name or "").lower()
    return any(t in n for t in _ONCOLOGY_TERMS)


# Anti-targets: proteins a drug commonly BINDS as a safety/ADME liability but never
# engages therapeutically. Their disease associations are toxicities (hERG->arrhythmia),
# not repurposing opportunities — so they must not GENERATE candidate indications. This
# only bites drugs whose target list came from bioactivity gap-fill (no curated
# mechanism); curated-mechanism drugs never list these. Extensible, deliberately small
# and high-confidence so we never drop a real therapeutic target.
_ANTI_TARGETS = frozenset({
    "KCNH2",   # hERG — cardiac repolarisation; the canonical cardiotoxicity anti-target
    "KCNQ1", "SCN5A",           # cardiac ion channels (QT/arrhythmia liability screens)
    "CYP3A4", "CYP2D6", "CYP2C9", "CYP1A2", "CYP2C19",  # drug-metabolising enzymes (ADME)
    "HTR2B",   # 5-HT2B — valvulopathy off-target liability
})


def _repurposing_targets(drug_genes: List[str]) -> List[str]:
    """The drug's targets minus safety anti-targets, for candidate GENERATION. Falls back
    to the full list if filtering would leave nothing (don't blank an anti-target-only drug)."""
    kept = [g for g in drug_genes if g.upper() not in _ANTI_TARGETS]
    return kept or drug_genes


def _is_oncology_drug(info: Dict) -> bool:
    """Is the QUERY drug itself an oncology drug — i.e. is its established use cancer?
    We look only at its top development-phase tier (its approved/lead indications, not
    the long tail of exploratory trials) and ask whether those are predominantly cancer.
    Palbociclib (phase-4 = breast cancer) -> True, so we repurpose it OUT of oncology;
    Nintedanib (top tier = IPF/scleroderma/NSCLC, ~half fibrosis) -> False, so oncology
    stays a valid repurposing target. Unknown/undocumented -> False (don't over-filter)."""
    inds = info.get("known_indications") or []
    if not inds:
        return False
    top = max((float(i.get("max_phase") or 0) for i in inds), default=0.0)
    if top <= 0:
        return False
    tier = [i for i in inds if float(i.get("max_phase") or 0) >= top - 0.001]
    onc = sum(1 for i in tier if _is_oncology_name(i.get("name", "")))
    return onc / max(1, len(tier)) >= 0.75


def screen_indications_for_drug(
    drug: str,
    area_filter: Optional[str] = None,
    max_candidates: int = MAX_SCORED_CANDIDATES,
    exclude_oncology: Optional[bool] = None,
) -> Dict:
    """
    Find and rank NEW candidate indications for a molecule.

    exclude_oncology: drop cancer indications. None (default) = DRUG-AWARE: exclude
    oncology only if the query drug is ITSELF an oncology drug (repurpose it outward);
    for a non-oncology drug (e.g. nintedanib, a fibrosis drug) oncology stays a valid
    repurposing target. Pass True/False to force. Ignored when area_filter is oncology.

    Candidates are ranked by a mechanism EVIDENCE score (target overlap, OT
    association, pathway, PPI), not by the drug's clinical/regulatory status
    (which is constant for a fixed drug). The 505(b)(2) feasibility is reported
    separately per candidate.
    """
    info = resolve_drug(drug)
    chembl_id = info["chembl_id"]
    drug_genes = info["targets"]

    if not chembl_id:
        return {"drug": drug, "error": "Could not resolve drug in ChEMBL", "candidates": []}

    # F7 biologic-aware path (audited 2026-07): classify modality once. For an antibody /
    # protein / peptide / oligo, the structure-based science stack (docking, PBPK, quantum,
    # CNS-MPO, med-chem liability, potency lead-viability) does NOT apply — mark it as an
    # explicit, first-class N/A rather than a silent null, and rank on the target/pathway
    # + trial/biomarker axis the biologic actually acts through.
    try:
        from services.modality import classify as _classify_modality
        _modality = _classify_modality(name=info.get("name", drug),
                                       smiles=info.get("smiles", ""),
                                       molecule_type=info.get("molecule_type", ""))
    except Exception as e:
        logger.debug(f"modality classify failed: {e}")
        _modality = {"modality": "unknown", "is_small_molecule": bool(info.get("smiles")),
                     "optimization": "smiles", "label": "unknown"}
    # A molecule is treated as small-molecule ONLY when we can actually analyse it as one
    # (a SMILES is present). An unknown-modality drug with NO SMILES (e.g. an antibody the
    # classifier could not label) must NOT be run through CNS-MPO / developability / potency
    # math — those would emit meaningless small-molecule numbers presented as real. Treat it
    # as non-small-molecule so those modules are shielded (marked N/A) instead.
    _has_smiles = bool((info.get("smiles") or "").strip())
    _is_small_molecule = _modality.get("is_small_molecule", True) and _has_smiles
    _is_biologic = not _is_small_molecule

    # Drug-aware oncology default: only an oncology drug repurposes OUT of oncology.
    onc_drug = _is_oncology_drug(info)
    if exclude_oncology is None:
        exclude_oncology = onc_drug

    cache_key = (f"drug:{drug.lower().strip()}|area:{(area_filter or '').lower()}"
                 f"|noonc:{int(exclude_oncology)}|v3")   # v3: + recovered-known-indication lane
    cache = _load_cache()
    if cache_key in cache:
        return cache[cache_key]

    if not drug_genes:
        return {"drug": info["name"], "chembl_id": chembl_id,
                "error": "No protein targets found for this molecule",
                "candidates": [], "drug_targets": []}

    # Repurposing target set = the drug's targets minus safety anti-targets (hERG etc.),
    # so a toxicity liability can't GENERATE or target-score a candidate indication.
    screen_genes = _repurposing_targets(drug_genes)
    excluded_anti = [g for g in drug_genes if g not in screen_genes]

    # 1. Candidate diseases from each target's Open Targets associations ────────
    # Fetched concurrently — each gene needs a live Ensembl lookup + an associations
    # query, and doing them serially for up to 12 targets was a large part of the
    # reverse screen's latency. The merge below stays serial (and order-stable).
    candidate_map: Dict[str, Dict] = {}

    def _gene_diseases(gene: str):
        ensembl = _gene_to_ensembl(gene)
        return (gene, _diseases_for_target(ensembl) if ensembl else [])

    with ThreadPoolExecutor(max_workers=8) as _ex:
        _gene_results = list(_ex.map(_gene_diseases, screen_genes[:12]))

    for gene, _diseases in _gene_results:
        for d in _diseases:
            if d["score"] < MIN_ASSOCIATION_SCORE:
                continue
            # Skip ontology SIGNS/PHENOTYPES (e.g. myelosuppression EFO_0007053,
            # HP_* terms) — these are toxicity consequences of hitting the target,
            # not diseases the drug could be repurposed INTO.
            if _is_phenotype_not_disease(d["efo_id"], d["therapeutic_areas"]):
                continue
            efo = d["efo_id"]
            if efo not in candidate_map:
                candidate_map[efo] = {
                    "disease":           d["name"],
                    "efo_id":            efo,
                    "therapeutic_areas": d["therapeutic_areas"],
                    "association_score": d["score"],
                    "via_target":        gene,
                    "sources":           {"genetic"},
                    "trial_count":       0,
                    "max_trial_phase":   0,
                    "lit_count":         0,
                }
            elif d["score"] > candidate_map[efo]["association_score"]:
                candidate_map[efo]["association_score"] = d["score"]
                candidate_map[efo]["via_target"] = gene

    # 1a-KG. KG-GENERATED candidates (gap #2 — the Rephetio/TxGNN inversion) ─────────
    # Instead of only scoring what OT target-associations surfaced, the validated DWPC
    # model scores the drug against EVERY Hetionet disease and adds its top P(treats)
    # indications — ones genetic association misses (imatinib→epilepsy / IgA-GN). The KG
    # GENERATES candidates (recall); the canonical composite still RANKS them (precision).
    # Bounded by Hetionet coverage (2016; sparse/newer drugs add nothing — fail-soft).
    try:
        from services.repurposing_predictor import top_diseases_for_drug
        _existing = {v["disease"].lower(): k for k, v in candidate_map.items()}
        for kg in top_diseases_for_drug(info.get("name", drug), k=25):
            nm = kg["disease"]; low = nm.lower()
            if low in _existing:                       # OT already has it → mark KG-corroborated
                v = candidate_map[_existing[low]]
                v["sources"] = set(v.get("sources", set())) | {"knowledge-graph"}
                v["kg_probability"] = kg["probability"]
                continue
            candidate_map[f"kg:{nm}"] = {
                "disease": nm, "efo_id": "", "therapeutic_areas": [],
                "association_score": 0.0, "via_target": "", "sources": {"knowledge-graph"},
                "trial_count": 0, "max_trial_phase": 0, "lit_count": 0,
                "kg_probability": kg["probability"],
            }
            _existing[low] = f"kg:{nm}"
    except Exception as e:
        logger.debug(f"KG candidate generation skipped: {e}")

    # 1a-DRKG. MODERN-KG candidates — same Rephetio/TxGNN inversion, but on DRKG
    # (~5.1k diseases / ~23k compounds, 6 sources) instead of 2016 Hetionet (213/8894).
    # This is what lets newer drugs Hetionet never saw get a KG signal — e.g. Palbociclib
    # (CDK4/6i, 2015) surfaces Rheumatoid Arthritis, which Hetionet cannot. Diseases carry
    # a MeSH/DOID id; the canonical composite still ranks them. Fail-soft on missing model.
    try:
        from services.drkg_predictor import top_diseases_for_drug as _drkg_top
        _existing = {v["disease"].lower(): k for k, v in candidate_map.items()}
        for kg in _drkg_top(info.get("name", drug), chembl_id=info.get("chembl_id", ""),
                            k=30, min_p=0.95):
            nm = kg["disease"]; low = nm.lower()
            if low in _existing:
                v = candidate_map[_existing[low]]
                v["sources"] = set(v.get("sources", set())) | {"knowledge-graph"}
                v["kg_probability"] = max(v.get("kg_probability", 0.0), kg["probability"])
                continue
            # Tag oncology from the disease name so exclude_oncology can filter DRKG
            # cancer candidates (they carry no OT therapeutic_area to key off otherwise).
            candidate_map[f"drkg:{nm}"] = {
                "disease": nm, "efo_id": "",
                "therapeutic_areas": (["cancer or benign tumor"] if _is_oncology_name(nm) else []),
                "association_score": 0.0, "via_target": "", "sources": {"knowledge-graph"},
                "trial_count": 0, "max_trial_phase": 0, "lit_count": 0,
                "kg_probability": kg["probability"], "mesh_id": kg.get("identifier", ""),
            }
            _existing[low] = f"drkg:{nm}"
    except Exception as e:
        logger.debug(f"DRKG candidate generation skipped: {e}")

    # 1a-PrimeKG. BROAD-COVERAGE candidates over PrimeKG's 17,080 diseases (the graph TxGNN
    # is built on) — far wider than Hetionet (213) or DRKG (~5.1k). Ranked by the validated
    # FUSION generator (DistMult recall pool -> treats-classifier + R-GCN, RRF; beats the
    # classifier-only cascade on the held-out harness). It is a candidate-GENERATOR (recall);
    # the canonical composite + disease-value + contraindication gate still RANK/filter them,
    # so tail noise sinks. A high min_score keeps only confident P(treats) leads. Fail-soft:
    # needs torch + pkg_gnn.pt (present in .venv312); degrades to the cascade otherwise.
    try:
        from services.primekg_predictor import top_diseases_for_drug as _pkg_top
        _existing = {v["disease"].lower(): k for k, v in candidate_map.items()}
        for kg in _pkg_top(info.get("name", drug), chembl_id=info.get("chembl_id", ""),
                           k=30, min_score=0.90):
            nm = kg["disease"]; low = nm.lower()
            if low in _existing:
                v = candidate_map[_existing[low]]
                v["sources"] = set(v.get("sources", set())) | {"knowledge-graph"}
                v["kg_probability"] = max(v.get("kg_probability", 0.0), kg["score"])
                continue
            candidate_map[f"pkg:{nm}"] = {
                "disease": nm, "efo_id": "",
                "therapeutic_areas": (["cancer or benign tumor"] if _is_oncology_name(nm) else []),
                "association_score": 0.0, "via_target": "", "sources": {"knowledge-graph"},
                "trial_count": 0, "max_trial_phase": 0, "lit_count": 0,
                "kg_probability": kg["score"],
            }
            _existing[low] = f"pkg:{nm}"
    except Exception as e:
        logger.debug(f"PrimeKG candidate generation skipped: {e}")

    # 1b. Add candidates from clinical trials + literature (real-world repurposing
    #     signal that genetic associations miss — e.g. an eye-drop trial of an
    #     oncology drug). Resolve each to a disease/EFO and merge by EFO id.
    def _merge_clinical(name: str, *, trial_count: int = 0,
                        max_trial_phase: int = 0, lit_count: int = 0,
                        outcome_signal: float = 0.0, failed_efficacy: int = 0,
                        failed_safety: int = 0, serious_ae: dict = None,
                        endpoint_verdict: str = None, endpoint_note: str = ""):
        # Safety trials can list a drug's symptom/AE as a studied "condition"
        # (e.g. a trial of nausea prophylaxis); the generic denylist keeps those
        # bare symptom terms out of the candidate indications.
        if _is_ae_symptom(name):
            return
        meta = _resolve_disease_meta(name)
        efo = meta.get("efo", "")
        if not efo:
            return
        # A toxicity term arriving via trials/literature can resolve to a
        # phenotype EFO (e.g. myelosuppression → EFO_0007053, areas ['phenotype']);
        # drop it here too so it can't seed a candidate through the clinical path.
        if _is_phenotype_not_disease(efo, meta.get("areas", [])):
            return
        e = candidate_map.get(efo)
        if e is None:
            e = candidate_map[efo] = {
                "disease": meta["name"], "efo_id": efo,
                "therapeutic_areas": meta["areas"], "association_score": 0.0,
                "via_target": "", "sources": set(),
                "trial_count": 0, "max_trial_phase": 0, "lit_count": 0,
                "trial_outcome_signal": 0.0, "failed_efficacy": 0, "failed_safety": 0,
            }
        if trial_count:
            e["sources"].add("clinical_trial")
            e["trial_count"]     = max(e["trial_count"], trial_count)
            e["max_trial_phase"] = max(e["max_trial_phase"], max_trial_phase)
            e["trial_outcome_signal"] = outcome_signal
            e["failed_efficacy"] = failed_efficacy
            e["failed_safety"]   = failed_safety
            if serious_ae is not None:
                e["serious_ae"] = serious_ae      # denominatored serious-AE for the Safety column
            # Display-only copy of the parsed primary-endpoint verdict for the Efficacy
            # column. Kept under a SEPARATE key so it never feeds _evidence_tier / the rank
            # (that path still reads c["endpoint_verdict"], which the merge does not set).
            if endpoint_verdict:
                e["ct_efficacy_verdict"] = endpoint_verdict
                e["ct_efficacy_note"] = endpoint_note or ""
        if lit_count:
            e["sources"].add("literature")
            e["lit_count"] = max(e["lit_count"], lit_count)

    # Stamp-at-fetch marker for the Verification stage: set once ClinicalTrials.gov
    # actually answers, so a per-candidate clinicaltrials `checked` flag is truthful.
    _ctgov_probe: Dict = {}
    trials: Dict[str, Dict] = {}
    try:
        trials = _trials_for_drug(info["name"], probe=_ctgov_probe)
        for v in sorted(trials.values(), key=lambda x: -x["count"])[:MAX_TRIAL_CONDITIONS]:
            _merge_clinical(v["name"], trial_count=v["count"], max_trial_phase=v["max_phase"],
                            outcome_signal=v.get("outcome_signal", 0.0),
                            failed_efficacy=v.get("failed_efficacy", 0),
                            failed_safety=v.get("failed_safety", 0),
                            serious_ae=v.get("serious_ae"),
                            endpoint_verdict=v.get("endpoint_verdict"),
                            endpoint_note=v.get("endpoint_note", ""))
    except Exception as e:
        logger.debug(f"trial merge failed: {e}")
    try:
        lit = _literature_for_drug(info["name"])
        for name, cnt in sorted(lit.items(), key=lambda x: -x[1])[:MAX_LIT_DISEASES]:
            _merge_clinical(name, lit_count=cnt)
    except Exception as e:
        logger.debug(f"literature merge failed: {e}")

    # 2. Exclude only APPROVED indications (the genuine novelty disqualifier).
    #    Prior development at a lower phase is kept — for a 505(b)(2)/hybrid play
    #    an already-de-risked indication is an asset, not a duplicate to hide; it
    #    surfaces as a candidate and is flagged as prior clinical art downstream.
    approved_inds = [k for k in info["known_indications"]
                     if float(k.get("max_phase") or 0) >= 4]
    known_names = [k["name"] for k in approved_inds]
    known_efos  = {k["efo_id"] for k in approved_inds if k["efo_id"]}

    def _is_approved_here(c: Dict) -> bool:
        # (1) Exact ontology-id match (fast path, when the ids happen to align).
        if c["efo_id"] and c["efo_id"] in known_efos:
            return True
        # (2) Order-independent, cross-ontology name match. ChEMBL stores approved
        #     indications as MeSH headings ("Arthritis, Rheumatoid", EFO id) while
        #     candidates arrive from Open Targets as canonical names ("rheumatoid
        #     arthritis", MONDO id) — so the id sets and the raw strings differ for
        #     the SAME disease. Token-set equivalence catches it (this is what let a
        #     drug's own approved indication, e.g. Baricitinib→RA, top the list).
        if any(_same_disease(c["disease"], kn) for kn in known_names):
            return True
        # (3) FDA-label supplement — catches indications ChEMBL still records below
        #     phase 4 but that are already FDA-approved (ChEMBL's per-indication
        #     phase lags real approvals). Label text is cached per drug (one fetch
        #     per screen). Fail-soft.
        try:
            from services.approval_supplement import is_label_approved_for
            return is_label_approved_for(info["name"], c["disease"])
        except Exception:
            return False

    # Final AE/symptom guard: a bare symptom / adverse-reaction term (Headache,
    # Nausea, Pruritus…) is never a repurposing INDICATION, no matter which
    # source proposed it — the genetic target→disease path can surface the same
    # symptom terms that the literature/trials paths do. Generic denylist, not a
    # per-drug list, so it stays consistent with the runtime-derived design.
    # F4 (audited 2026-07): narrow BROAD umbrella candidates (vasculitis) to the drug-
    # relevant subtype (EGPA for an IL-5 drug) BEFORE the approval filter — so a card names
    # the EXACT disease, and an approved subtype is then correctly excluded. Bounded to the
    # strongest few (each subtype resolve is a live lookup); memoized.
    _narrowed = 0
    for c in candidate_map.values():
        if _narrowed >= 8:
            break
        if float(c.get("association_score") or 0) < 0.20:
            continue
        nb = narrow_broad_disease(c.get("disease", ""), drug_genes)
        if nb and nb.lower() != (c.get("disease", "") or "").lower():
            c["narrowed_from"] = c.get("disease")
            c["disease"] = nb
            c["efo_id"] = ""
            _narrowed += 1
    candidates = [c for c in candidate_map.values()
                  if not _is_approved_here(c) and not _is_ae_symptom(c["disease"])]

    # RECOVERED KNOWN INDICATIONS (mechanism validation, not a novelty claim) ─────
    # An approved indication is normally excluded from the NOVEL candidate list above.
    # But when the mechanism machinery INDEPENDENTLY GENERATED it — target→disease
    # association / knowledge-graph / a real target link, i.e. it arrived through the
    # same evidence path a novel lead would — that is a validation signal: the platform
    # re-derived a real approved use from mechanism alone. Silently dropping it hides
    # exactly the proof the ranking works. Surface it in a SEPARATE, clearly-labelled
    # lane (never mixed into the novel leads, never presented as a repurposing claim).
    # Fully general: the test is "was it mechanistically generated?", with no per-drug
    # or per-disease logic — any drug whose on-label use the engine recovers qualifies.
    def _mechanistically_generated(c: Dict) -> bool:
        srcs = c.get("sources", set())
        return (float(c.get("association_score") or 0) > 0
                or "genetic" in srcs or "knowledge-graph" in srcs
                or bool(c.get("via_target")))
    recovered_candidates = [c for c in candidate_map.values()
                            if _is_approved_here(c)
                            and not _is_ae_symptom(c["disease"])
                            and _mechanistically_generated(c)]

    # 3a. Exclude oncology unless the user explicitly asked for it ──────────────
    if exclude_oncology and not _area_is_oncology(area_filter or ""):
        candidates = [c for c in candidates if not _is_oncology_area(c["therapeutic_areas"])]

    # 3b. Optional therapeutic-area filter (data-driven from Open Targets) ──────
    if area_filter:
        candidates = [c for c in candidates if _matches_area(c["therapeutic_areas"], area_filter)]

    # 4. Tractability: downrank ultra-rare genetic leaves (e.g. neonatal syndromes),
    #    then drop diseases with no known drugs — undruggable / too rare to develop.
    def _is_genetic_leaf(c: Dict) -> bool:
        return any("genetic, familial or congenital" in ta.lower() for ta in c["therapeutic_areas"])

    def _rank_score(c: Dict) -> float:
        # Credit real-world signal so trial/literature leads (genetic assoc = 0)
        # are not trimmed away before scoring.
        base = c["association_score"]
        srcs = c.get("sources", set())
        if "clinical_trial" in srcs:
            base = max(base, 0.5 + 0.1 * c.get("max_trial_phase", 0))
        elif "literature" in srcs:
            base = max(base, 0.25)
        # KG-generated candidates (genetic assoc = 0) rank on the validated DWPC
        # P(treats) so they survive trimming and get scored alongside OT candidates.
        if "knowledge-graph" in srcs:
            base = max(base, 0.4 * c.get("kg_probability", 0.0))
        return base * (0.5 if _is_genetic_leaf(c) else 1.0)

    candidates.sort(key=_rank_score, reverse=True)
    pool = candidates[: max_candidates + 20]
    # Tractability lookups run concurrently — one live Open Targets call per candidate,
    # previously serial over ~45 candidates (a dominant contributor to the >5 min hangs).
    def _known(c):
        return _disease_known_drugs(c["efo_id"]) if c.get("efo_id") else 0
    with ThreadPoolExecutor(max_workers=8) as _ex:
        _kd = list(_ex.map(_known, pool))
    tractable: List[Dict] = []
    for c, kd in zip(pool, _kd):
        c["known_drugs"] = kd
        # keep druggable indications AND strong KG-generated leads (no efo to look up)
        if c["known_drugs"] > 0 or c.get("kg_probability", 0.0) >= 0.5:
            tractable.append(c)
    # If the tractability filter is too aggressive, fall back to the ranked pool
    candidates = (tractable or pool)
    candidates.sort(key=_rank_score, reverse=True)
    candidates = candidates[:max_candidates]

    drug_compound = {
        "chembl_id":   chembl_id,
        "name":        info["name"],
        "max_phase":   info["max_phase"],
        "indications": "; ".join(k["name"] for k in info["known_indications"]),
        "targets":     ";".join(drug_genes),
    }

    # CNS-MPO (Wager 2010): can this molecule cross the blood-brain barrier? This
    # gates whether a CNS candidate indication is even physically plausible.
    from services.cns_mpo import cns_mpo as _cns_mpo
    # CNS-MPO is a small-molecule (MW/TPSA/HBD/logP) property score — not applicable to a
    # biologic. Mark it N/A explicitly instead of returning a null-filled shell.
    if _is_biologic:
        drug_mpo = {"applicable": False, "cns_druggable": None,
                    "note": f"Not applicable — {_modality.get('label', 'biologic')} (CNS-MPO scores small-molecule properties)."}
    else:
        drug_mpo = _cns_mpo(props=info.get("properties") or {}, smiles=info.get("smiles") or "")

    def _is_cns_disease(areas, name):
        # Classify from the candidate's REAL Open Targets therapeutic areas (no
        # keyword guessing), so psychiatric/CNS diseases are never mislabelled.
        try:
            from services.therapeutic_context import therapeutic_context
            return therapeutic_context(name, therapeutic_areas=areas).get("is_cns", False)
        except Exception:
            blob = (" ".join(areas or []) + " " + (name or "")).lower()
            return any(h in blob for h in ("nervous system", "central nervous", "brain",
                                           "neuro", "psychiat", "mental", "cognit"))

    try:
        from services.repurposing_engine import score_compound_for_disease
        from services.disease_ontology import resolve_disease, get_ppi_network
    except Exception as e:
        return {"drug": info["name"], "chembl_id": chembl_id,
                "error": f"Scoring dependencies unavailable: {e}", "candidates": []}

    try:
        from services import developability as _dev
    except Exception:
        _dev = None

    # Direction + appropriateness inputs, computed ONCE per drug ────────────────
    from services.disease_appropriateness import infer_drug_action, appropriateness
    # A. Drug action for the direction engine: curated ChEMBL mechanism, else inferred
    #    from the name (nintedanib has no curated mechanism, so direction was blind).
    try:
        from services.repurposing_engine import _actions_for_molecules
        _curated_actions = _actions_for_molecules([chembl_id]).get(chembl_id, {}) or {}
    except Exception:
        _curated_actions = {}
    _inferred_action = infer_drug_action(info.get("name", drug), screen_genes)
    drug_action_map = dict(_curated_actions)
    if not drug_action_map and _inferred_action:
        drug_action_map = {g: _inferred_action for g in screen_genes}
    drug_action = _inferred_action
    if not drug_action and _curated_actions:
        _blob = " ".join(_curated_actions.values()).lower()
        drug_action = ("inhibitor" if any(w in _blob for w in ("inhib", "antagon", "block"))
                       else "agonist" if "agonist" in _blob else "")
    # F7: a monoclonal antibody / antisense oligo NEUTRALIZES its target by default
    # (mepolizumab blocks IL-5, omalizumab blocks IgE), so the direction is antagonist
    # even without a curated ChEMBL action — otherwise direction/appropriateness runs blind.
    if not drug_action and _is_biologic and _modality.get("modality") in ("antibody", "oligonucleotide"):
        drug_action = "inhibitor"
        if not drug_action_map:
            drug_action_map = {g: "inhibitor" for g in screen_genes}
    # C. The drug's serious FAERS reactions — to flag AE-not-target candidates.
    # Stamp-at-fetch marker for the Verification stage (drug-level, shared by candidates).
    _faers_probe: Dict = {}
    try:
        from services.safety_filter import _faers_serious_reactions
        _faers_rx = _faers_serious_reactions(info.get("name", drug), probe=_faers_probe)
    except Exception:
        _faers_rx = []
    _faers_total = sum(cnt for _, cnt in _faers_rx)

    drug_gene_set = {g.upper() for g in drug_genes}
    # per-indication development phase (drug's phase in THIS disease) → phase-scale prior
    # clinical art in the scorer, so a merely-STUDIED indication is not credited like an
    # APPROVED one (the confounding that ranked tried/failed indications at the top).
    _known_inds = info.get("known_indications", [])
    # APPROVED indications (Phase 4) — drives Layer-2 off-label cannibalization: a candidate
    # in the same organ as one of these is already prescribable off-label.
    _approved_inds = [k.get("name", "") for k in _known_inds if float(k.get("max_phase") or 0) >= 4]
    # Layer 4 — base-compound exclusivity profile (drug-level, computed once): Orange Book
    # patents (small molecule) or Purple Book cliff (biologic).
    try:
        from services.regulatory_505b2 import exclusivity_profile as _excl_profile
        _excl_505b2 = _excl_profile(info.get("name", drug), _is_biologic)
    except Exception as _ee:
        logger.debug(f"exclusivity profile skipped: {_ee}")
        _excl_505b2 = {}
    def _ind_phase_for(dis: str):
        return indication_phase_for(dis, _known_inds)
    # Each candidate is scored independently (its own DB reads + live trial/patent/lit
    # enrichment), so the candidates are scored concurrently on a bounded pool. This is the
    # dominant cost of a cold reverse screen. Worker count is kept at/under the DB connection
    # pools' headroom (raised to 8) so a pool is never exhausted into silent-empty results.
    def _score_one(c):
        _attach_trials_by_name(c, trials)   # so a failed trial reaches a mechanism-only candidate
        dinfo = resolve_disease(c["disease"]) or {}
        disease_genes    = [t["gene_symbol"] for t in dinfo.get("targets", [])[:40]]
        disease_pathways = dinfo.get("pathways", [])
        ppi_adj = get_ppi_network(disease_genes[:15]) if disease_genes else {}

        disease_weights = {t["gene_symbol"].upper():
                           (t.get("quality_score") or t.get("genetic_score") or t.get("score", 0.0))
                           for t in dinfo.get("targets", []) if t.get("gene_symbol")}
        # FIX 1 — pass the SAME therapeutic_areas + studied_for_disease the canonical/dossier
        # path now passes, so the reverse LIST card and the dossier score the pair identically
        # (single source of truth; otherwise the canonical edit would open a new list-vs-dossier
        # drift on the LoF/congenital gate and the FAERS confounding guard).
        _studied_here = False
        try:
            from services.repurposing_scorer import approved_chembls_for_disease
            if chembl_id:
                _studied_here = chembl_id in approved_chembls_for_disease(
                    [chembl_id], [], c["disease"], min_phase=1)
        except Exception:
            _studied_here = False
        sr = score_compound_for_disease(
            drug_compound, c["disease"], disease_genes, disease_pathways, ppi_adj, screen_genes,
            disease_gene_weights=disease_weights or None, drug_actions=drug_action_map or None,
            trial_count=c.get("trial_count", 0), max_trial_phase=c.get("max_trial_phase", 0),
            trial_outcome=c.get("trial_outcome_signal", 0.0),
            indication_phase=_ind_phase_for(c["disease"]),
            therapeutic_areas=c.get("therapeutic_areas") or [],
            studied_for_disease=_studied_here,
        )
        sc = sr["scores"]

        # Disease-appropriateness gate (B loss-of-function/developmental, C adverse-event):
        # would the drug WORSEN or does it CAUSE this disease? Gates ranking + actionability;
        # does not overwrite the mechanistic composite.
        appr = appropriateness(info.get("name", drug), c["disease"], c["therapeutic_areas"],
                               drug_action, faers_reactions=_faers_rx, faers_total=_faers_total)
        overlap = sorted(drug_gene_set & {g.upper() for g in disease_genes})

        # Clinical/literature signal SPECIFIC to this indication (varies per disease,
        # unlike the drug's global phase): existing trials of this drug here + papers.
        trial_count = c.get("trial_count", 0)
        lit_count   = c.get("lit_count", 0)
        clinical_signal = min(1.0,
            0.5 * (c.get("max_trial_phase", 0) / 4.0)
            + 0.3 * min(1.0, trial_count / 5.0)
            + 0.2 * min(1.0, lit_count / 15.0)
        )

        evidence = (
            EVIDENCE_WEIGHTS["target"]      * sc.get("target", 0.0)
            + EVIDENCE_WEIGHTS["association"] * c["association_score"]
            + EVIDENCE_WEIGHTS["pathway"]     * sc.get("pathway", 0.0)
            + EVIDENCE_WEIGHTS["ppi"]         * sc.get("ppi", 0.0)
            + EVIDENCE_WEIGHTS["clinical"]    * clinical_signal
        )

        # Plain-language "why this might treat it"
        bits = []
        if overlap:
            bits.append(
                f"{info['name'].title()} acts on {', '.join(overlap[:3])}, which Open Targets "
                f"associates with {c['disease']} (association {c['association_score']:.2f})."
            )
        elif c["association_score"] > 0:
            bits.append(
                f"Open Targets links {info['name'].title()}'s target {c['via_target']} to "
                f"{c['disease']} (association {c['association_score']:.2f})."
            )
        if trial_count:
            ph = c.get("max_trial_phase", 0)
            bits.append(
                f"{trial_count} clinical trial{'s' if trial_count > 1 else ''} of {info['name'].title()} "
                f"in {c['disease']}{f' (up to Phase {ph})' if ph else ''} already exist."
            )
        if lit_count and not trial_count:
            bits.append(f"Co-mentioned with {info['name'].title()} in {lit_count} PubMed records.")
        if not bits:
            bits.append(f"Weak, association-only lead for {c['disease']}.")
        is_cns = _is_cns_disease(c["therapeutic_areas"], c["disease"])
        if is_cns and drug_mpo.get("score") is not None:
            bits.append(
                f"BBB risk: CNS-MPO {drug_mpo['score']}/6 - brain penetration uncertain "
                f"for this CNS indication." if drug_mpo["score"] < 4.0
                else f"Brain-penetrant (CNS-MPO {drug_mpo['score']}/6).")
        rationale = " ".join(bits)

        # Area-aware developability: judge the molecule on the properties that
        # matter for THIS indication's route (skin / ocular / CNS / oral).
        developability = {}
        if _dev is not None and info.get("smiles"):
            try:
                developability = _dev.score(
                    info["smiles"], area=area_filter or "",
                    therapeutic_areas=c["therapeutic_areas"],
                )
            except Exception:
                pass

        # Adaptive actionable cutoff (applied AFTER approved-indication exclusion, so
        # a trial link can never re-admit an on-label use). Base bar is relaxed when
        # there is real-world trial movement in THIS indication (human evidence
        # overrides raw computational score) or when the disease is orphan (rare
        # diseases carry sparse data, so a strict bar blinds us to real leads).
        # CTPA Rule 1 (cohesion) + Rule 2 (registry ghost) are now applied inside the
        # canonical scorer (score_compound_for_disease), so every surface inherits
        # them; here we just read the verdicts for the card chips.
        ctpa = sr.get("ctpa", {"phantom": False, "multiplier": 1.0})
        ctpa_registry = sr.get("registry", {"ghost": False, "multiplier": 1.0})
        # PrimeKG ground-truth cross-check: the x0.15 kill multiplier is already folded into
        # the composite by the canonical scorer; here we read the verdict so a labeled
        # contraindication can force the candidate non-actionable and be explained in why_not.
        pkg = sr.get("primekg", {}) or {}
        pkg_contra = pkg.get("relation") == "contraindication"

        try:
            _orphan = bool(is_orphan_candidate(c["disease"], c.get("efo_id", "")))
        except Exception:
            _orphan = False
        _has_trial = trial_count > 0 or c.get("max_trial_phase", 0) > 0

        # Actionability is judged on the RAW composite score against the well-separated
        # 2026-07 bands (base 0.40 = the "Promising" floor). The calibrated percentile is
        # reported as background context (see the "calibration" field) but is NOT the gate:
        # the signal/noise gap widened so far that the null collapsed near zero and the
        # percentile SATURATES (raw 0.15 and 0.83 both read top ~98%), so gating on it made
        # nearly every mechanistically-linked candidate "actionable". score_calibration.py
        # abandoned the percentile for tiering for exactly this reason; we match it here.
        # Adaptive relaxations: a trial link in this indication and an orphan disease each
        # lower the raw bar.
        eff_cutoff = _BASE_ACTIONABLE_CUTOFF
        if _has_trial:
            eff_cutoff = min(eff_cutoff, _TRIAL_LINK_CUTOFF)
        if _orphan:
            eff_cutoff *= _ORPHAN_CUTOFF_FACTOR
        eff_cutoff = round(eff_cutoff, 3)
        actionable = sr["composite_score"] >= eff_cutoff
        cutoff_kind = "raw"

        entry = {
            **c,
            "sources":              sorted(c.get("sources", set())),   # set → JSON-safe list
            "evidence_score":       round(evidence, 4),
            "score":                sr["composite_score"],   # canonical score (same everywhere)
            "composite_score":      sr["composite_score"],
            "safety":               sr.get("safety", {}),     # toxicity cross-filter verdict
            "coverage":             sr.get("coverage", {}),    # target-coverage gate verdict
            "clinical_constraints": sr.get("clinical_constraints", {}),  # CCH severity/chronicity fit
            "trial_failure":        sr.get("trial_failure", {}),          # P3 demonstrated efficacy/safety failure
            "directional_evidence": sr.get("directional_evidence", {}),   # P3 typed drug->gene->disease triples
            "ctpa":                 ctpa,                                 # CTPA Rule 1 cohesion gate verdict
            "registry":             ctpa_registry,                        # CTPA Rule 2 ghost audit verdict
            "mechanism_scope":      sr.get("mechanism_scope", {}),  # target-mediated? (else score N/A)
            "calibration":          sr.get("calibration", {}),      # percentile/tier vs null background
            "primekg":              pkg,   # ground-truth relation + direction-aware P(treats) cross-check
            "plausibility":         _plausibility_for(info.get("name", drug), c["disease"]),  # validated DWPC P(treats), where covered
            "lead_viability":       _lead_viability_for(c.get("chembl_id", ""), info.get("name", drug),
                                                        c["disease"], overlap[:10] or drug_genes,
                                                        smiles=info.get("smiles", "")),  # potency funnel (window deferred to dossier)
            "effective_cutoff":     eff_cutoff,
            "cutoff_kind":          cutoff_kind,
            "actionable":           bool(actionable and appr["appropriate"] and not pkg_contra),
            "appropriateness":      appr,   # B/C: would the drug worsen or cause this disease?
            "evidence_tier":        _evidence_tier(c),   # real-world strength: supported / promising / contradicted
            # Human-readable Safety / Efficacy columns sourced from ClinicalTrials.gov
            # (serious-AE rate with denominator; primary-endpoint verdict). DISPLAY ONLY —
            # not folded into composite_score or any rank key (no double counting). A pair
            # with no trials reads "No trial data" for BOTH, never a safe/effective default.
            "safety_ct":            _ct_safety_col(c.get("serious_ae"),
                                                   n_trials=trial_count,
                                                   phase=c.get("max_trial_phase", 0)),
            "efficacy_ct":          _ct_efficacy_col(c.get("ct_efficacy_verdict") or c.get("endpoint_verdict"),
                                                     n_trials=trial_count,
                                                     phase=c.get("max_trial_phase", 0),
                                                     note=c.get("ct_efficacy_note") or c.get("endpoint_note", "")),
            "scores":               sc,
            "clinical_signal":      round(clinical_signal, 4),
            "trial_count":          trial_count,
            "max_trial_phase":      c.get("max_trial_phase", 0),
            "lit_count":            lit_count,
            "overlapping_targets":  overlap[:10],
            "disease_gene_count":   len(disease_genes),
            "rationale":            rationale,
            "developability":       developability,
            "developability_score": developability.get("score"),
            "cns_relevant":         is_cns,
            "bbb_caution":          bool(is_cns and drug_mpo.get("score") is not None
                                         and drug_mpo["score"] < 4.0),
            "exploratory":          c.get("known_drugs", 0) == 0,
            # Would a pharma company pursue THIS disease? (burden × unmet-need × market)
            "disease_value":        _disease_value_for(c["disease"], c.get("efo_id", "")),
        }
        # Layer 2 — commercial & off-label friction (bounded ranking adjustment + flags)
        try:
            from services.commercial_friction import assess as _comm_assess
            entry["commercial_friction"] = _comm_assess(
                c["disease"], _approved_inds, disease_value=entry.get("disease_value"))
        except Exception as _ce:
            logger.debug(f"commercial friction skipped: {_ce}")
            entry["commercial_friction"] = {}
        # Layer 4 — 505(b)(2) feasibility & exclusivity (route + runway + formulation whitespace)
        try:
            from services.regulatory_505b2 import feasibility as _reg_feas
            entry["regulatory_505b2"] = _reg_feas(
                c["disease"], info.get("name", drug), _excl_505b2, _is_biologic,
                smiles=info.get("smiles", ""), disease_value=entry.get("disease_value"))
        except Exception as _re:
            logger.debug(f"505b2 feasibility skipped: {_re}")
            entry["regulatory_505b2"] = {}

        # ── Coverage-aware Verification stage (label / lane only) ─────────────────
        # Record which negative sources were actually consulted for THIS pair. Uses the
        # drug-level stamp-at-fetch markers (_ctgov_probe / _faers_probe) so a genuinely
        # queried-but-clean source is never confused with one that was never run. Writes
        # nothing to the composite and adds no rank multiplier: it only relabels a
        # mechanism-only tier to a visibly distinct "unverified" lane at the SAME rank
        # factor, so a sparse novel candidate is flagged, never sunk below junk.
        try:
            from services import verification as _verif
            _ipf = _ind_phase_for(c["disease"])
            _approved = bool(_ipf is not None and int(_ipf) >= 4)
            entry["evidence_balance"] = _verif.evidence_balance(sr, {
                "faers_queried":      bool(_faers_probe.get("faers_queried")),
                "faers_suppressed":   _approved,
                "faers_total":        _faers_total,
                "primekg_suppressed": _approved,
                "coverage_suppressed": _approved,
                "ctgov_queried":      bool(_ctgov_probe.get("ctgov_queried")),
                "trial_count":        int(c.get("trial_count", 0) or 0),
                "endpoint_verdict":   c.get("endpoint_verdict"),
            })
            if entry["evidence_balance"].get("verdict") == "insufficient_coverage":
                _et = entry.get("evidence_tier") or {}
                if _et.get("tier") in ("preclinical", "mechanistic"):
                    entry["evidence_tier"] = {
                        "tier": "unverified", "label": "Unverified hypothesis",
                        "note": ("Negative sources not yet checked for this pair; "
                                 "hypothesis only.")}
        except Exception as _ve:
            logger.debug(f"verification stage skipped: {_ve}")
            entry["evidence_balance"] = {}

        entry["why_not"] = _why_not(entry)
        return entry

    with ThreadPoolExecutor(max_workers=4) as _score_pool:
        scored: List[Dict] = list(_score_pool.map(_score_one, candidates))

    # DEFAULT SORT = pharma repurposing value × pair evidence. The disease Value Score
    # (burden/unmet-need/market-fit) modulates the mechanistic evidence so that a great
    # drug-fit for a TRIVIAL indication (headache, skin allergy) sinks below a good fit
    # for a serious, high-unmet-need one — the fix for "that's not what a pharma company
    # would repurpose for". Unmapped diseases get a neutral 0.5 (never zeroed on a data
    # gap); association strength remains a separate, visible evidence field.
    # Real-world evidence tier nudges rank order: a trial-supported candidate outranks an
    # equal-scoring mechanism-only one, and a contradicted one sinks — without overwriting
    # the mechanistic composite (which stays the displayed score).
    # FIX 2 — pair-specific + DIRECTIONAL evidence factor. ABSENCE of clinical/trial art is
    # NEUTRAL (exactly 1.0), never a penalty — a genuinely novel mechanism-strong indication must
    # not be buried behind an equal-mechanism already-tried one purely for prior art. Positive
    # human evidence keeps a >1.0 edge but only as a MILD TIE-BREAKER (<=5%), not a dominating
    # multiplier that reorders mechanism. A FAILED / contradicted trial FOR THIS indication stays
    # negative (<1.0). lit_count (PubMed co-mention) no longer carries any ranking weight: the
    # "literature" tier now reads 1.0 (its co-mention signal remains a displayed field only, and
    # the reverse `evidence_score` axis that also carried it was already display-only — not in _rank).
    _TIER_RANK = {"trial-supported": 1.05, "tested-unverified": 1.02, "promising": 1.01,
                  # absent / mechanism-only / preclinical / literature / coverage-unverified:
                  # NEUTRAL — novelty is relabelled to the hypothesis lane, never rank-demoted.
                  "literature": 1.00, "preclinical": 1.00, "mechanistic": 1.00,
                  "unverified": 1.00,
                  # negative-directional evidence FOR THIS indication stays below 1.0:
                  "biomarker-only": 0.80, "failed-endpoint": 0.60, "contradicted": 0.55}
    def _rank(x):
        dv = x.get("disease_value") or {}
        vw = dv.get("value_score", 0.5) if dv else 0.5
        appr_f = (x.get("appropriateness") or {}).get("factor", 1.0)
        tier_f = _TIER_RANK.get((x.get("evidence_tier") or {}).get("tier", ""), 1.0)
        comm_f = float((x.get("commercial_friction") or {}).get("multiplier", 1.0))  # Layer 2
        return float(x.get("composite_score", 0.0)) * (0.4 + 0.6 * vw) * appr_f * tier_f * comm_f
    scored.sort(key=_rank, reverse=True)

    # F8 (audited 2026-07): collapse duplicate candidates that resolve to the SAME disease
    # under different EFO ids / name variants (e.g. two 'pulmonary hypertension, primary'
    # rows from the OT and trial sources). The sort above put the best-ranked instance
    # first, so keep first-seen per normalized disease key and drop the rest.
    _seen_dis, _dedup = set(), []
    for _c in scored:
        _dn = (_c.get("disease") or "").lower()
        _key = " ".join(sorted(w for w in re.split(r"[^a-z0-9]+", _dn) if len(w) > 3)) or _dn
        if _key in _seen_dis:
            continue
        _seen_dis.add(_key)
        _dedup.append(_c)
    scored = _dedup

    # Score the RECOVERED known indications through the IDENTICAL scorer, so their
    # mechanistic evidence chain is computed the same way a novel lead's is. They are
    # kept in a separate list (never merged into `scored`) and each is flagged as an
    # approved use recovered by mechanism, so the UI can present it as validation, not
    # as a new repurposing opportunity.
    recovered_scored: List[Dict] = []
    if recovered_candidates:
        with ThreadPoolExecutor(max_workers=4) as _rec_pool:
            recovered_scored = list(_rec_pool.map(_score_one, recovered_candidates))
        for _e in recovered_scored:
            _e["is_known_indication"]    = True
            _e["recovered_by_mechanism"] = True
            _e["novelty"]                = "known indication, recovered by mechanism"
        recovered_scored.sort(key=_rank, reverse=True)
        _seen_rec, _rec_dedup = set(), []
        for _c in recovered_scored:
            _dn = (_c.get("disease") or "").lower()
            _key = " ".join(sorted(w for w in re.split(r"[^a-z0-9]+", _dn) if len(w) > 3)) or _dn
            if _key in _seen_rec:
                continue
            _seen_rec.add(_key)
            _rec_dedup.append(_c)
        recovered_scored = _rec_dedup

    # Positive-Pivot generation: for CCH-crushed candidates, turn the mismatch into a
    # discovery lead (severe/orphan variant, dose-sparing combination, or brain-
    # penetrant analogue) instead of only eliminating it. Bounded to the strongest
    # few crushed candidates (each pivot does a live OT/ChEMBL lookup). Hypotheses.
    try:
        from services.positive_pivots import generate_pivots
        _pivots_done = 0
        for e in scored:
            if _pivots_done >= 5:
                break
            cc = e.get("clinical_constraints") or {}
            if not cc.get("penalized"):
                continue
            e["pivots"] = generate_pivots(
                info["name"], e["disease"], cch=cc, smiles=info.get("smiles", ""),
                disease_genes=[t for t in e.get("overlapping_targets", [])] or drug_genes)
            _pivots_done += 1
    except Exception as e:
        logger.debug(f"pivot generation skipped: {e}")

    # Flag candidates that are a KNOWN but sub-phase-4 (studied, not approved) indication —
    # prior clinical art kept intentionally, NOT a novel discovery. This is why e.g.
    # eosinophilic esophagitis (Phase 2 for mepolizumab) legitimately appears as a
    # candidate; the card labels it so, and it is NOT listed as "excluded".
    _studied = [(k["name"], float(k.get("max_phase") or 0)) for k in info["known_indications"]
                if 0 < float(k.get("max_phase") or 0) < 4]
    for _c in scored:
        for _nm, _ph in _studied:
            if _same_disease(_c.get("disease", ""), _nm):
                _c["prior_studied_phase"] = _ph
                break

    result = {
        "drug":              info["name"],
        "chembl_id":         chembl_id,
        "smiles":            info.get("smiles", ""),
        "drug_targets":      drug_genes[:20],
        # Only APPROVED indications are actually excluded from results; studied (sub-
        # phase-4) ones remain candidates (flagged prior_studied_phase), so split them.
        "excluded_indications": [k["name"] for k in info["known_indications"]
                                 if float(k.get("max_phase") or 0) >= 4],
        "studied_indications":  [k["name"] for k in info["known_indications"]
                                 if 0 < float(k.get("max_phase") or 0) < 4],
        "known_indications": [k["name"] for k in info["known_indications"]],
        "area_filter":       area_filter or "",
        "exclude_oncology":  exclude_oncology,
        "is_oncology_drug":  onc_drug,
        "excluded_anti_targets": excluded_anti,
        "candidate_count":   len(scored),
        "candidates":        scored,
        # Approved indications the mechanism engine independently RE-DERIVED (validation
        # signal, not novel leads). Scored by the same composite; ranked, labelled.
        "recovered_indications": recovered_scored,
        "cns_mpo":           drug_mpo,
        # F7 biologic-aware: expose modality + an explicit list of the small-molecule-only
        # modules that do NOT apply, so the UI shows first-class "N/A (biologic)" instead
        # of blank/zero panels a pharma reviewer would read as a gap.
        # Provenance + freshness: the data lineage behind these candidates, each source
        # dated + integrity-scored so a reviewer sees exactly where the ranking came from
        # (and that e.g. the ChEMBL bioactivity snapshot is aging while trials are live).
        "provenance":        _reverse_provenance(),
        "modality":          _modality,
        "not_applicable":    ([
            "Molecular docking (no small-molecule structure)",
            "PBPK exposure (small-molecule perfusion model)",
            "Quantum chemistry (small-molecule electronic structure)",
            "CNS-MPO (small-molecule properties)",
            "Structure-aware med-chem liability",
            "Potency lead-viability (IC50/Ki funnel)",
        ] if _is_biologic else []),
        "modality_basis":    (f"{_modality.get('label', 'biologic')} — ranked on target/pathway "
                              f"and clinical/biomarker evidence; small-molecule structure modules "
                              f"are not applicable." if _is_biologic else ""),
        "_ts":               time.time(),
    }
    cache[cache_key] = result
    _save_cache(cache)
    return result
