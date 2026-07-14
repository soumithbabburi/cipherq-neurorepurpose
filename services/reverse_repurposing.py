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
from pathlib import Path
from typing import Dict, List, Optional

import requests  # noqa: F401

from services import http_client

logger = logging.getLogger(__name__)

OT_URL      = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CT_BASE     = "https://clinicaltrials.gov/api/v2/studies"
NCBI_BASE   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
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
    def is_orphan_candidate(_name: str) -> bool:     # fail-soft stub
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


def _gene_to_ensembl(symbol: str) -> str:
    """Resolve a gene symbol → Ensembl target ID via Open Targets search."""
    q = """
    query($q: String!) {
      search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
        hits { id name entity }
      }
    }"""
    hits = _ot_graphql(q, {"q": symbol}).get("search", {}).get("hits", [])
    return hits[0]["id"] if hits else ""


def _disease_known_drugs(efo_id: str) -> int:
    """Count of drugs + clinical candidates for a disease — a tractability signal."""
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


def _trials_for_drug(drug_name: str, page_size: int = 50) -> Dict[str, Dict]:
    """Conditions studied in ClinicalTrials.gov for this drug → candidate indications,
    WITH per-condition trial OUTCOMES (P3): completed / with-results / failed-for-
    efficacy / failed-for-safety / ongoing, and a directional outcome_signal
    (-1..+1). A drug that FAILED a trial for efficacy in an indication is real
    negative evidence, not merely "a trial exists".

    Relationship-aware: `query.intr` so the drug is an INTERVENTION, not merely
    mentioned (the v2 param is `query.intr`, not `query.intervention`).
    """
    out: Dict[str, Dict] = {}
    try:
        r = http_client.get(CT_BASE, params={"query.intr": drug_name,
                                          "pageSize": page_size, "format": "json"}, timeout=15)
        if r and r.ok:
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
                                           "ongoing": 0})
                    e["count"] += 1
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
            for e in out.values():
                e["outcome_signal"] = _trial_outcome_signal(e)
    except Exception as ex:
        logger.debug(f"trials fetch failed: {ex}")
    return out


def _trial_outcome_signal(e: Dict) -> float:
    """Directional trial-outcome signal in [-1, 1]: positive for completed / with-
    results trials, negative for efficacy/safety failures in this indication."""
    pos = 0.20 * e.get("completed", 0) + 0.35 * e.get("with_results", 0) + 0.08 * e.get("ongoing", 0)
    neg = 0.60 * e.get("failed_efficacy", 0) + 0.40 * e.get("failed_safety", 0)
    return round(max(-1.0, min(1.0, pos - neg)), 3)


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


def _diseases_for_target(ensembl_id: str, size: int = 250) -> List[Dict]:
    """Open Targets diseases associated with a target, with therapeutic areas."""
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


def resolve_drug(drug: str) -> Dict:
    """
    Resolve a drug name or ChEMBL ID to {chembl_id, name, max_phase, targets, known_indications}.
    Salts are mapped to their parent molecule. Targets and indications are fetched live from ChEMBL.
    """
    chembl_id = ""

    if re.fullmatch(r"CHEMBL\d+", drug.strip(), re.I):
        chembl_id = drug.strip().upper()
    else:
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

    # ChEMBL records mechanisms / indications inconsistently across the parent and
    # its salt forms — gather every form and union the results.
    parent_id = _to_parent(chembl_id)
    ids = _molecule_forms(parent_id)
    if chembl_id not in ids:
        ids.insert(0, chembl_id)

    # Display name / phase / SMILES: prefer the parent's record, fall back to the salt's.
    name, max_phase, smiles = drug, 0, ""
    for cid in (parent_id, chembl_id):
        d = _molecule_detail(cid)
        smiles = smiles or (d.get("molecule_structures") or {}).get("canonical_smiles", "")
        if d.get("pref_name"):
            name = d["pref_name"]
            max_phase = d.get("max_phase") or 0
            break

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
    known_indications: List[Dict] = list(by_label.values())

    return {
        "chembl_id":         parent_id,
        "name":              name,
        "max_phase":         max_phase,
        "smiles":            smiles,
        "targets":           targets,
        "known_indications": known_indications,
    }


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
    cache_key = f"pair:{(chembl_id or drug_name).lower()}|{disease.lower().strip()}{_pk}"
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
    sr = score_compound_for_disease(compound, disease, disease_genes, disease_pathways, ppi,
                                    drug_genes, drug_actions=drug_actions,
                                    disease_gene_weights=disease_weights or None,
                                    mechanistic_prior=mechanistic_prior)
    out = {"score": sr["composite_score"], "composite_score": sr["composite_score"],
           "scores": sr["scores"], "safety": sr.get("safety", {}),
           "coverage": sr.get("coverage", {}),
           "clinical_constraints": sr.get("clinical_constraints", {}),
           "trial_failure": sr.get("trial_failure", {}),
           "directional_evidence": sr.get("directional_evidence", {}),
           "calibration": sr.get("calibration", {}),
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
    # Accept a few common synonyms without hardcoding disease lists
    synonyms = {
        "dermatology": ["skin", "dermat", "integument"],
        "ophthalmology": ["eye", "ophthalm", "visual", "ocular"],
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


def screen_indications_for_drug(
    drug: str,
    area_filter: Optional[str] = None,
    max_candidates: int = MAX_SCORED_CANDIDATES,
    exclude_oncology: bool = True,
) -> Dict:
    """
    Find and rank NEW candidate indications for a molecule.

    exclude_oncology: drop cancer indications (default True) — repurposing an
    oncology drug into more oncology isn't a repurposing story. Ignored when the
    requested area_filter is itself oncology.

    Candidates are ranked by a mechanism EVIDENCE score (target overlap, OT
    association, pathway, PPI), not by the drug's clinical/regulatory status
    (which is constant for a fixed drug). The 505(b)(2) feasibility is reported
    separately per candidate.
    """
    cache_key = (f"drug:{drug.lower().strip()}|area:{(area_filter or '').lower()}"
                 f"|noonc:{int(exclude_oncology)}")
    cache = _load_cache()
    if cache_key in cache:
        return cache[cache_key]

    info = resolve_drug(drug)
    chembl_id = info["chembl_id"]
    drug_genes = info["targets"]

    if not chembl_id:
        return {"drug": drug, "error": "Could not resolve drug in ChEMBL", "candidates": []}
    if not drug_genes:
        return {"drug": info["name"], "chembl_id": chembl_id,
                "error": "No protein targets found for this molecule",
                "candidates": [], "drug_targets": []}

    # 1. Candidate diseases from each target's Open Targets associations ────────
    candidate_map: Dict[str, Dict] = {}
    for gene in drug_genes[:12]:
        ensembl = _gene_to_ensembl(gene)
        if not ensembl:
            continue
        for d in _diseases_for_target(ensembl):
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

    # 1b. Add candidates from clinical trials + literature (real-world repurposing
    #     signal that genetic associations miss — e.g. an eye-drop trial of an
    #     oncology drug). Resolve each to a disease/EFO and merge by EFO id.
    def _merge_clinical(name: str, *, trial_count: int = 0,
                        max_trial_phase: int = 0, lit_count: int = 0,
                        outcome_signal: float = 0.0, failed_efficacy: int = 0,
                        failed_safety: int = 0):
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
        if lit_count:
            e["sources"].add("literature")
            e["lit_count"] = max(e["lit_count"], lit_count)

    try:
        trials = _trials_for_drug(info["name"])
        for v in sorted(trials.values(), key=lambda x: -x["count"])[:MAX_TRIAL_CONDITIONS]:
            _merge_clinical(v["name"], trial_count=v["count"], max_trial_phase=v["max_phase"],
                            outcome_signal=v.get("outcome_signal", 0.0),
                            failed_efficacy=v.get("failed_efficacy", 0),
                            failed_safety=v.get("failed_safety", 0))
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
    candidates = [c for c in candidate_map.values()
                  if not _is_approved_here(c) and not _is_ae_symptom(c["disease"])]

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
        return base * (0.5 if _is_genetic_leaf(c) else 1.0)

    candidates.sort(key=_rank_score, reverse=True)
    pool = candidates[: max_candidates + 15]
    tractable: List[Dict] = []
    for c in pool:
        c["known_drugs"] = _disease_known_drugs(c["efo_id"])
        if c["known_drugs"] > 0:
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

    drug_gene_set = {g.upper() for g in drug_genes}
    scored: List[Dict] = []
    for c in candidates:
        dinfo = resolve_disease(c["disease"]) or {}
        disease_genes    = [t["gene_symbol"] for t in dinfo.get("targets", [])[:40]]
        disease_pathways = dinfo.get("pathways", [])
        ppi_adj = get_ppi_network(disease_genes[:15]) if disease_genes else {}

        disease_weights = {t["gene_symbol"].upper():
                           (t.get("quality_score") or t.get("genetic_score") or t.get("score", 0.0))
                           for t in dinfo.get("targets", []) if t.get("gene_symbol")}
        sr = score_compound_for_disease(
            drug_compound, c["disease"], disease_genes, disease_pathways, ppi_adj, drug_genes,
            disease_gene_weights=disease_weights or None,
            trial_count=c.get("trial_count", 0), max_trial_phase=c.get("max_trial_phase", 0),
            trial_outcome=c.get("trial_outcome_signal", 0.0),
        )
        sc = sr["scores"]
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

        try:
            _orphan = bool(is_orphan_candidate(c["disease"]))
        except Exception:
            _orphan = False
        _has_trial = trial_count > 0 or c.get("max_trial_phase", 0) > 0

        # Actionability is judged on the CALIBRATED percentile when the null exists
        # (a top-percentile lead is actionable even if its raw score is modest —
        # this reconciles the 0.40 raw bar with "top 1%"). The adaptive exceptions
        # relax the PERCENTILE bar: a trial link and an orphan disease each lower it.
        cal = sr.get("calibration", {})
        if cal.get("basis") == "calibrated" and cal.get("percentile") is not None:
            pct_bar = 0.85                      # top 15% = actionable
            if _has_trial:
                pct_bar = 0.75                  # trial link → top 25%
            if _orphan:
                pct_bar -= 0.10                 # orphan → relax further
            actionable = cal["percentile"] >= pct_bar
            eff_cutoff = round(pct_bar, 3)      # a percentile bar, shown as such
            cutoff_kind = "percentile"
        else:                                   # heuristic fallback (no null yet)
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
            "plausibility":         _plausibility_for(info.get("name", drug), c["disease"]),  # validated DWPC P(treats), where covered
            "lead_viability":       _lead_viability_for(c.get("chembl_id", ""), info.get("name", drug),
                                                        c["disease"], overlap[:10] or drug_genes,
                                                        smiles=info.get("smiles", "")),  # potency funnel (window deferred to dossier)
            "effective_cutoff":     eff_cutoff,
            "cutoff_kind":          cutoff_kind,
            "actionable":           bool(actionable),
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
        scored.append(entry)

    # DEFAULT SORT = pharma repurposing value × pair evidence. The disease Value Score
    # (burden/unmet-need/market-fit) modulates the mechanistic evidence so that a great
    # drug-fit for a TRIVIAL indication (headache, skin allergy) sinks below a good fit
    # for a serious, high-unmet-need one — the fix for "that's not what a pharma company
    # would repurpose for". Unmapped diseases get a neutral 0.5 (never zeroed on a data
    # gap); association strength remains a separate, visible evidence field.
    def _rank(x):
        dv = x.get("disease_value") or {}
        vw = dv.get("value_score", 0.5) if dv else 0.5
        return float(x.get("composite_score", 0.0)) * (0.4 + 0.6 * vw)
    scored.sort(key=_rank, reverse=True)

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

    result = {
        "drug":              info["name"],
        "chembl_id":         chembl_id,
        "smiles":            info.get("smiles", ""),
        "drug_targets":      drug_genes[:20],
        "known_indications": [k["name"] for k in info["known_indications"]],
        "area_filter":       area_filter or "",
        "exclude_oncology":  exclude_oncology,
        "candidate_count":   len(scored),
        "candidates":        scored,
        "cns_mpo":           drug_mpo,
        "_ts":               time.time(),
    }
    cache[cache_key] = result
    _save_cache(cache)
    return result
