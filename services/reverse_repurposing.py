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
import re
import time
from pathlib import Path
from typing import Dict, List, Optional

import requests

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


# ── Cache ───────────────────────────────────────────────────────────────────

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
        CACHE_FILE.write_text(json.dumps(cache), encoding="utf-8")
    except Exception:
        pass


# ── Open Targets helpers ────────────────────────────────────────────────────

def _ot_graphql(query: str, variables: dict, timeout: int = 12) -> dict:
    try:
        r = requests.post(OT_URL, json={"query": query, "variables": variables},
                          timeout=timeout, headers={"Content-Type": "application/json"})
        if r.ok:
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


def _trials_for_drug(drug_name: str, page_size: int = 50) -> Dict[str, Dict]:
    """Conditions studied in ClinicalTrials.gov for this drug → candidate indications."""
    out: Dict[str, Dict] = {}
    try:
        r = requests.get(CT_BASE, params={"query.term": drug_name,
                                          "pageSize": page_size, "format": "json"}, timeout=15)
        if r.ok:
            for s in r.json().get("studies", []):
                ps = s.get("protocolSection", {})
                conds = ps.get("conditionsModule", {}).get("conditions", []) or []
                phases = ps.get("designModule", {}).get("phases", []) or []
                ph = _max_phase_from_phases(phases)
                for cond in conds:
                    k = cond.strip().lower()
                    if not k:
                        continue
                    e = out.setdefault(k, {"name": cond.strip(), "count": 0, "max_phase": 0})
                    e["count"] += 1
                    e["max_phase"] = max(e["max_phase"], ph)
    except Exception as e:
        logger.debug(f"trials fetch failed: {e}")
    return out


def _literature_for_drug(drug_name: str, retmax: int = 120) -> Dict[str, int]:
    """Disease MeSH headings co-occurring with the drug in PubMed → candidate indications."""
    counts: Dict[str, int] = {}
    try:
        sr = requests.get(f"{NCBI_BASE}/esearch.fcgi",
                          params={"db": "pubmed", "term": f'"{drug_name}"[tiab]',
                                  "retmax": retmax, "retmode": "json"}, timeout=10)
        ids = sr.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return counts
        fr = requests.get(f"{NCBI_BASE}/efetch.fcgi",
                          params={"db": "pubmed", "id": ",".join(ids), "retmode": "xml"},
                          timeout=25)
        if fr.ok:
            import xml.etree.ElementTree as ET
            root = ET.fromstring(fr.content)
            for mh in root.iter("MeshHeading"):
                d = mh.find("DescriptorName")
                if d is None or not (d.text or "").strip():
                    continue
                name = d.text.strip()
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
        r = requests.get(f"{CHEMBL_BASE}/molecule/{chembl_id}.json", timeout=10)
        if r.ok:
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
        r = requests.get(f"{CHEMBL_BASE}/molecule_form.json",
                         params={"parent_chembl_id": parent_id, "limit": 50, "format": "json"},
                         timeout=10)
        if r.ok:
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
            r = requests.get(f"{CHEMBL_BASE}/molecule/search.json",
                             params={"q": drug, "limit": 25}, timeout=15)
            if r.ok:
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

    known_indications: List[Dict] = []
    seen = set()
    for cid in ids:
        try:
            r = requests.get(f"{CHEMBL_BASE}/drug_indication.json",
                             params={"molecule_chembl_id": cid, "limit": 100,
                                     "format": "json"}, timeout=10)
            if r.ok:
                for ind in r.json().get("drug_indications", []):
                    label = ind.get("mesh_heading") or ind.get("efo_term") or ""
                    efo   = (ind.get("efo_id") or "").replace(":", "_")
                    if label and label.lower() not in seen:
                        seen.add(label.lower())
                        known_indications.append({"name": label, "efo_id": efo})
        except Exception:
            pass

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
                         drug_name: str = "") -> Dict:
    """THE repurposing score for a (drug, disease) pair — used identically by the
    Repurpose card, the analysis header, and the Evidence Dossier so the number is
    the same everywhere. Composite of the 6 screens; cached per pair (deterministic)."""
    cache_key = f"pair:{(chembl_id or drug_name).lower()}|{disease.lower().strip()}"
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
    sr = score_compound_for_disease(compound, disease, disease_genes, disease_pathways, ppi, drug_genes)
    out = {"score": sr["composite_score"], "composite_score": sr["composite_score"],
           "scores": sr["scores"], "drug_genes": drug_genes[:20], "_ts": time.time()}
    cache[cache_key] = out
    _save_cache(cache)
    return out


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
                        max_trial_phase: int = 0, lit_count: int = 0):
        meta = _resolve_disease_meta(name)
        efo = meta.get("efo", "")
        if not efo:
            return
        e = candidate_map.get(efo)
        if e is None:
            e = candidate_map[efo] = {
                "disease": meta["name"], "efo_id": efo,
                "therapeutic_areas": meta["areas"], "association_score": 0.0,
                "via_target": "", "sources": set(),
                "trial_count": 0, "max_trial_phase": 0, "lit_count": 0,
            }
        if trial_count:
            e["sources"].add("clinical_trial")
            e["trial_count"]     = max(e["trial_count"], trial_count)
            e["max_trial_phase"] = max(e["max_trial_phase"], max_trial_phase)
        if lit_count:
            e["sources"].add("literature")
            e["lit_count"] = max(e["lit_count"], lit_count)

    try:
        trials = _trials_for_drug(info["name"])
        for v in sorted(trials.values(), key=lambda x: -x["count"])[:MAX_TRIAL_CONDITIONS]:
            _merge_clinical(v["name"], trial_count=v["count"], max_trial_phase=v["max_phase"])
    except Exception as e:
        logger.debug(f"trial merge failed: {e}")
    try:
        lit = _literature_for_drug(info["name"])
        for name, cnt in sorted(lit.items(), key=lambda x: -x[1])[:MAX_LIT_DISEASES]:
            _merge_clinical(name, lit_count=cnt)
    except Exception as e:
        logger.debug(f"literature merge failed: {e}")

    # 2. Exclude the drug's existing indications (only surface NEW ones) ────────
    known_names = {k["name"].lower() for k in info["known_indications"]}
    known_efos  = {k["efo_id"] for k in info["known_indications"] if k["efo_id"]}

    def _is_known(c: Dict) -> bool:
        if c["efo_id"] in known_efos:
            return True
        cn = c["disease"].lower()
        return any(cn == kn or cn in kn or kn in cn for kn in known_names if len(kn) > 4)

    candidates = [c for c in candidate_map.values() if not _is_known(c)]

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

        sr = score_compound_for_disease(
            drug_compound, c["disease"], disease_genes, disease_pathways, ppi_adj, drug_genes,
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

        entry = {
            **c,
            "sources":              sorted(c.get("sources", set())),   # set → JSON-safe list
            "evidence_score":       round(evidence, 4),
            "score":                sr["composite_score"],   # canonical score (same everywhere)
            "composite_score":      sr["composite_score"],
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
            "exploratory":          c.get("known_drugs", 0) == 0,
        }
        scored.append(entry)

    # Tractable indications (those with existing drugs to reference for 505(b)(2)) rank
    # above untreatable rare leaves; within each tier, by the canonical repurposing score.
    scored.sort(key=lambda x: (x.get("known_drugs", 0) > 0, x["composite_score"]), reverse=True)

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
        "_ts":               time.time(),
    }
    cache[cache_key] = result
    _save_cache(cache)
    return result
