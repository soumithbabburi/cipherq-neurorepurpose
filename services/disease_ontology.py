"""
Disease Ontology Service
Resolves diseases → gene targets, pathways, MeSH IDs using:
  - Open Targets Platform GraphQL (primary, no key needed)
  - NCBI eUtils (MeSH resolution)
  - STRING DB (protein-protein interactions)
  - Reactome (pathway data via Open Targets)
All results cached in data/ontology_cache.json.
"""
import json, logging, time
from pathlib import Path
import requests

logger = logging.getLogger(__name__)

CACHE_FILE = Path(__file__).parent.parent / "data" / "ontology_cache.json"
OT_URL     = "https://api.platform.opentargets.org/api/v4/graphql"
NCBI_URL   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
STRING_URL = "https://string-db.org/api/json"

# Known neurological/CNS disease prevalence (US patients)
# Used for orphan-drug threshold (< 200,000 = eligible)
PREVALENCE = {
    "alzheimer": 5_800_000, "parkinson": 1_000_000,
    "multiple sclerosis": 1_000_000, "ms": 1_000_000,
    "epilepsy": 3_200_000, "schizophrenia": 2_000_000,
    "bipolar": 5_000_000, "depression": 21_000_000,
    "als": 30_000, "amyotrophic lateral sclerosis": 30_000,
    "huntington": 30_000, "huntington disease": 30_000,
    "huntingtons": 30_000,
    "narcolepsy": 200_000, "tuberous sclerosis": 50_000,
    "rett": 15_000, "rett syndrome": 15_000,
    "friedreich ataxia": 15_000, "spinal muscular atrophy": 15_000,
    "duchenne": 20_000, "dravet": 15_000,
    "niemann pick": 10_000, "gaucher": 20_000,
    "fabry": 10_000, "prion": 5_000,
    "creutzfeldt jakob": 450,
}

# MONDO/EFO ID hints for common diseases (verified against Open Targets Platform)
EFO_HINTS = {
    "alzheimer":                "MONDO_0004975",
    "parkinson":                "MONDO_0005180",
    "als":                      "MONDO_0004976",
    "amyotrophic lateral sclerosis": "MONDO_0004976",
    "multiple sclerosis":       "MONDO_0005301",
    "ms":                       "MONDO_0005301",
    "huntington":               "MONDO_0007739",
    "epilepsy":                 "HP_0001250",
    "schizophrenia":            "MONDO_0005090",
    "depression":               "MONDO_0002050",
    "bipolar":                  "MONDO_0004985",
    "autism":                   "MONDO_0005260",
    "adhd":                     "MONDO_0007977",
    "anxiety":                  "MONDO_0011762",
    "stroke":                   "MONDO_0005098",
    "migraine":                 "MONDO_0005277",
}


def _load_cache() -> dict:
    try:
        if CACHE_FILE.exists():
            with open(CACHE_FILE, encoding="utf-8") as f:
                return json.load(f)
    except Exception:
        pass
    return {}


def _save_cache(cache: dict):
    try:
        CACHE_FILE.parent.mkdir(exist_ok=True)
        with open(CACHE_FILE, "w", encoding="utf-8") as f:
            json.dump(cache, f)
    except Exception:
        pass


def _ot_graphql(query: str, variables: dict, timeout: int = 12) -> dict:
    try:
        r = requests.post(OT_URL, json={"query": query, "variables": variables},
                          timeout=timeout, headers={"Content-Type": "application/json"})
        if r.ok:
            return r.json().get("data", {})
    except Exception as e:
        logger.debug(f"Open Targets GraphQL error: {e}")
    return {}


def resolve_disease(disease_name: str) -> dict:
    """
    Resolve a disease name to Open Targets disease record.
    Returns: {disease_id, name, description, targets: [{gene_symbol, score}], pathways: []}
    Cached per disease name.
    """
    key = disease_name.lower().strip()
    cache = _load_cache()
    cache_key = f"resolve:{key}"
    if cache_key in cache:
        return cache[cache_key]

    # Try MONDO/EFO hint first (fast path for common diseases)
    efo_id = None
    for hint_key, hint_id in EFO_HINTS.items():
        if hint_key in key or key in hint_key:
            efo_id = hint_id
            break

    # Search Open Targets if no hint or hint returned nothing
    def _search_ot(term: str):
        search_q = """
        query($q: String!) {
          search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 3}) {
            hits { id name entity }
          }
        }"""
        data = _ot_graphql(search_q, {"q": term})
        hits = data.get("search", {}).get("hits", [])
        return hits[0]["id"] if hits else None

    if not efo_id:
        efo_id = _search_ot(disease_name)

    if not efo_id:
        return {}

    # Get disease info + top targets + pathways
    detail_q = """
    query($id: String!) {
      disease(efoId: $id) {
        id name description
        associatedTargets(page: {index: 0, size: 40}) {
          rows {
            target {
              approvedSymbol approvedName id
              pathways { pathway pathwayId }
            }
            score
          }
        }
      }
    }"""
    data = _ot_graphql(detail_q, {"id": efo_id})
    disease = data.get("disease", {})
    if not disease:
        # Hint ID may be stale — try search fallback
        fallback_id = _search_ot(disease_name)
        if fallback_id and fallback_id != efo_id:
            data = _ot_graphql(detail_q, {"id": fallback_id})
            disease = data.get("disease", {})
            if disease:
                efo_id = fallback_id
    if not disease:
        return {}

    rows = disease.get("associatedTargets", {}).get("rows", [])

    targets = [
        {
            "gene_symbol": r["target"]["approvedSymbol"],
            "gene_name":   r["target"]["approvedName"],
            "ensembl_id":  r["target"]["id"],
            "score":       round(r["score"], 4),
        }
        for r in rows
        if r.get("score", 0) > 0.05
    ]

    # Build pathway map
    pathway_map = {}
    for r in rows:
        if r.get("score", 0) < 0.1:
            continue
        for pw in r["target"].get("pathways", []):
            pid = pw.get("pathwayId", "")
            if not pid:
                continue
            if pid not in pathway_map:
                pathway_map[pid] = {"pathway_id": pid,
                                    "pathway_name": pw.get("pathway", ""),
                                    "genes": []}
            pathway_map[pid]["genes"].append(r["target"]["approvedSymbol"])

    pathways = sorted(pathway_map.values(),
                      key=lambda x: len(x["genes"]), reverse=True)[:20]

    result = {
        "disease_id":   disease["id"],
        "disease_name": disease.get("name", disease_name),
        "description":  (disease.get("description") or "")[:300],
        "targets":      targets[:40],
        "pathways":     pathways,
    }
    cache[cache_key] = result
    _save_cache(cache)
    return result


def get_disease_gene_set(disease_name: str, top_n: int = 30) -> list:
    """Return top gene symbols for disease, sorted by evidence score."""
    info = resolve_disease(disease_name)
    return [t["gene_symbol"] for t in
            sorted(info.get("targets", []), key=lambda x: x.get("score", 0), reverse=True)[:top_n]]


def get_disease_pathways(disease_name: str) -> list:
    """Return Reactome pathways associated with disease."""
    return resolve_disease(disease_name).get("pathways", [])


def is_orphan_candidate(disease_name: str) -> bool:
    """True if disease prevalence < 200,000 US patients (FDA orphan threshold)."""
    name = disease_name.lower()
    for k, v in PREVALENCE.items():
        if k in name or name in k:
            return v < 200_000
    return False


def get_ppi_network(gene_symbols: list, species: int = 9606) -> dict:
    """
    Fetch protein-protein interaction network from STRING DB.
    Returns {gene: [interacting_genes]} adjacency.
    """
    if not gene_symbols:
        return {}
    cache = _load_cache()
    cache_key = f"ppi:{','.join(sorted(gene_symbols[:10]))}"
    if cache_key in cache:
        return cache[cache_key]

    try:
        r = requests.get(f"{STRING_URL}/network",
            params={"identifiers": "%0d".join(gene_symbols[:20]),
                    "species": species, "required_score": 700,
                    "caller_identity": "neurorepurpose_platform"},
            timeout=10)
        if not r.ok:
            return {}

        adj = {}
        for link in r.json():
            a, b = link.get("preferredName_A",""), link.get("preferredName_B","")
            if a and b:
                adj.setdefault(a, []).append(b)
                adj.setdefault(b, []).append(a)

        cache[cache_key] = adj
        _save_cache(cache)
        return adj
    except Exception as e:
        logger.debug(f"STRING PPI error: {e}")
        return {}


def resolve_mesh_id(disease_name: str) -> str:
    """Resolve disease name to MeSH ID via NCBI eUtils."""
    cache = _load_cache()
    cache_key = f"mesh:{disease_name.lower().strip()}"
    if cache_key in cache:
        return cache[cache_key]
    try:
        r = requests.get(f"{NCBI_URL}/esearch.fcgi",
            params={"db": "mesh", "term": disease_name, "retmax": 1, "retmode": "json"},
            timeout=8)
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if ids:
            r2 = requests.get(f"{NCBI_URL}/esummary.fcgi",
                params={"db": "mesh", "id": ids[0], "retmode": "json"}, timeout=8)
            docs = r2.json().get("result", {})
            uid = ids[0]
            mesh_id = docs.get(uid, {}).get("ds_meshui", "")
            if mesh_id:
                cache[cache_key] = mesh_id
                _save_cache(cache)
                return mesh_id
    except Exception as e:
        logger.debug(f"MeSH resolve error: {e}")
    return ""
