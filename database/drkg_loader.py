"""
DRKG Neurological Disease Subset Loader
Dynamically discovers disease entities from the graph itself,
resolves their names, filters to neurological diseases by keyword,
then builds a compact cache for the knowledge graph visualiser.

Source: https://github.com/gnn4dr/DRKG
Usage:  python -m database.drkg_loader
        python -m database.drkg_loader C:\\path\\to\\drkg.tar.gz
Output: data/drkg_neuro.json
"""

import gzip
import io
import json
import logging
import re
import tarfile
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import requests

logger = logging.getLogger(__name__)

DRKG_URL    = "https://dgl-data.s3-us-west-2.amazonaws.com/dataset/DRKG/drkg.tar.gz"
CACHE_FILE  = Path(__file__).parent.parent / "data" / "drkg_neuro.json"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
NCBI_BASE   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OLS_BASE    = "https://www.ebi.ac.uk/ols/api"

# ── Neurological keyword filter ────────────────────────────────────────────────
# Any disease whose resolved name contains one of these (case-insensitive) is kept
NEURO_KEYWORDS = {
    "alzheimer", "parkinson", "multiple sclerosis", "amyotrophic",
    "huntington", "epilep", "schizophr", "depress", "bipolar",
    "autis", "autism", "migraine", "attention deficit", "adhd",
    "neurodegenera", "dementia", "lewy body", "frontotemporal",
    "motor neuron", "spinal muscular", "spinocerebellar",
    "narcolepsy", "tourette", "rett ", "dravet", "prion",
    "creutzfeldt", "friedreich", "niemann", "gaucher",
    "traumatic brain", "stroke", "cerebrovascular",
    "glioblas", "glioma", "brain tumor", "meningi",
    "neuropath", "peripheral neuropath",
}

# ── Relation type patterns ─────────────────────────────────────────────────────
_TREAT_KEYS  = {"CtD", "CDtDO", "indication", "Pa", "T", "treat"}
_DG_KEYS     = {"DaG", "DuG", "DdG", "DrG", "disease_gene"}
_TARGET_KEYS = {
    "CbG", "CdG", "CuG", "CiG",
    "INHIBITOR", "ACTIVATOR", "AGONIST", "ANTAGONIST",
    "BINDER", "BLOCKER", "MODULATOR", "SUBSTRATE",
    "target", "enzyme", "transporter", "carrier",
    "pharmacologically_interacts",
}
_PPI_KEYS    = {"BINDING", "REACTION", "CATALYSIS", "INTERACTS", "PHYSICAL"}


def _rel_type(rel: str) -> str:
    """Classify DRKG relation. Format: 'Source::Code::HeadType:TailType'"""
    parts = rel.split("::")
    # Extract the relation code (middle component) — do exact matches to avoid
    # false positives from letter "T" appearing in "HETIONET", "DISEASE", etc.
    rc = parts[1].upper() if len(parts) >= 2 else rel.upper()
    full = rel.upper()

    # Treat (Compound → Disease)
    if rc in {"T", "CTD", "PA", "CDTDO", "TREATS", "PALLIATES", "INDICATION"}:
        return "treat"
    if any(k in full for k in ("TREATS", "INDICATION", "PALLIATES")):
        return "treat"

    # Disease → Gene
    if rc in {"DAG", "DUG", "DDG", "DRG", "G+", "G-", "MD", "J"}:
        return "dis_gene"
    if "DISEASE_GENE" in full:
        return "dis_gene"

    # Compound → Gene targets
    if rc in {"CBG", "CDG", "CUG", "CIG", "CCG"}:
        return "target"
    if any(k in full for k in (
            "INHIBITOR", "ACTIVATOR", "AGONIST", "ANTAGONIST",
            "BINDER", "BLOCKER", "MODULATOR", "SUBSTRATE",
            "PHARMACOLOGICALLY", "ENZYME", "TRANSPORTER", "CARRIER")):
        return "target"

    # PPI
    if rc in {"BINDING", "REACTION", "CATALYSIS", "INTERACTS", "PHYSICAL"}:
        return "ppi"
    if any(k in full for k in ("BINDING", "REACTION", "CATALYSIS")):
        return "ppi"

    return "other"


def _etype(entity: str) -> str:
    return entity.split("::")[0] if "::" in entity else "Unknown"


def _esrc(entity: str) -> str:
    parts = entity.split("::")
    if len(parts) > 1:
        # Handles both "DrugBank" and "MESH:D009369" → returns "MESH"
        return parts[1].split(":")[0]
    return ""


def _eid(entity: str) -> str:
    parts = entity.split("::")
    if len(parts) >= 3:
        return "::".join(parts[2:])
    elif len(parts) == 2:
        # e.g. "MESH:D009369" → "D009369"
        sub = parts[1].split(":", 1)
        return sub[1] if len(sub) > 1 else sub[0]
    return entity


# ── File loading ───────────────────────────────────────────────────────────────

def _decode_content(raw: bytes) -> List[str]:
    if raw[:2] == b"\x1f\x8b":          # inner gzip
        raw = gzip.decompress(raw)
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    return text.splitlines()


def _extract_from_tar(fileobj) -> List[str]:
    with tarfile.open(fileobj=fileobj, mode="r:*") as tf:
        members = tf.getmembers()
        real = [m for m in members
                if not m.name.startswith("._") and m.size > 1000]
        member = (
            next((m for m in real if m.name.endswith("drkg.tsv")), None)
            or next((m for m in real if m.name.endswith(".tsv")),    None)
            or next((m for m in real if m.name.endswith(".tsv.gz")), None)
        )
        if not member:
            raise RuntimeError(
                f"No .tsv in archive. Members: {[m.name for m in members]}")
        print(f"  Member: {member.name}  ({member.size // (1024*1024)} MB)")
        raw = tf.extractfile(member).read()
    return _decode_content(raw)


def _load_tsv_from_file(path: Path) -> List[str]:
    suffix = "".join(path.suffixes).lower()
    print(f"  Local file: {path}  ({path.stat().st_size // (1024*1024)} MB)")
    if suffix in (".tar.gz", ".tgz"):
        lines = _extract_from_tar(io.BytesIO(path.read_bytes()))
    else:
        lines = _decode_content(path.read_bytes())
    print(f"  Loaded {len(lines):,} triplets.")
    return lines


def _download_tsv() -> List[str]:
    print("Downloading DRKG from S3…")
    resp = requests.get(DRKG_URL, stream=True, timeout=300)
    resp.raise_for_status()
    chunks, total = [], 0
    for chunk in resp.iter_content(chunk_size=512 * 1024):
        chunks.append(chunk)
        total += len(chunk)
        mb = total // (1024 * 1024)
        if mb % 10 == 0 and mb > 0:
            print(f"  … {mb} MB")
    print(f"  Done ({total // (1024*1024)} MB). Extracting…")
    lines = _extract_from_tar(io.BytesIO(b"".join(chunks)))
    print(f"  Extracted {len(lines):,} triplets.")
    return lines


# ── Phase 1: discover all disease entities that appear in edges ────────────────

def _discover_disease_entities(lines: List[str]) -> Tuple[
    Counter,            # disease entity → edge count
    List[dict],         # raw treat  edges  {h, r, t}
    List[dict],         # raw dg     edges
    List[dict],         # raw target edges
    Set[str],           # all compound entity IDs seen
]:
    """Single pass over all lines — collect everything we need."""
    dis_count:   Counter   = Counter()
    all_cmp_ids: Set[str]  = set()
    treat_all:   List[dict] = []
    dg_all:      List[dict] = []
    target_all:  List[dict] = []

    total     = len(lines)
    milestone = max(total // 10, 1)

    for i, line in enumerate(lines):
        if i % milestone == 0:
            print(f"  {i * 100 // total}%…", end=" ", flush=True)

        parts = line.strip().split("\t")
        if len(parts) != 3:
            continue
        h, r, t = parts
        ht, tt = _etype(h), _etype(t)
        rtype  = _rel_type(r)

        if rtype == "treat":
            if ht == "Compound" and tt == "Disease":
                treat_all.append({"h": h, "r": r, "t": t})
                dis_count[t] += 1
                all_cmp_ids.add(h)
            elif ht == "Disease" and tt == "Compound":
                treat_all.append({"h": t, "r": r, "t": h})
                dis_count[h] += 1
                all_cmp_ids.add(t)

        elif rtype == "dis_gene":
            if ht == "Disease" and tt == "Gene":
                dg_all.append({"h": h, "r": r, "t": t})
                dis_count[h] += 1
            elif ht == "Gene" and tt == "Disease":
                dg_all.append({"h": t, "r": r, "t": h})
                dis_count[t] += 1

        elif rtype == "target":
            # Collect all — filter to neuro compounds AFTER Phase 1 completes,
            # so we don't miss edges that appear before their compound's treat edge
            if ht == "Compound" and tt == "Gene":
                target_all.append({"h": h, "r": r, "t": t})
            elif ht == "Gene" and tt == "Compound":
                target_all.append({"h": t, "r": r, "t": h})

    print()   # newline after progress dots
    return dis_count, treat_all, dg_all, target_all, all_cmp_ids


# ── Phase 2: resolve disease names ────────────────────────────────────────────

def _resolve_one_disease(entity_id: str) -> Tuple[str, str]:
    """Return (entity_id, human-readable name). Uses appropriate API per ID type."""
    src = _esrc(entity_id).lower()
    eid = _eid(entity_id)

    # ── Name is embedded directly in entity string ─────────────────────────────
    # e.g. Disease::Disease::Alzheimer Disease
    #      Disease::BioPortal::Alzheimer_disease
    if src in ("disease", "bioportal", "name", "label"):
        return entity_id, eid.replace("_", " ")

    # ── DOID → EBI OLS API ─────────────────────────────────────────────────────
    if "doid" in src:
        doid_num = eid.lstrip("0") or "0"
        obo_id   = f"DOID:{eid}"
        try:
            r = requests.get(f"{OLS_BASE}/terms",
                params={"obo_id": obo_id}, timeout=8)
            if r.ok:
                terms = r.json().get("_embedded", {}).get("terms", [])
                if terms:
                    return entity_id, terms[0].get("label", eid)
        except Exception:
            pass
        return entity_id, f"DOID:{eid}"

    # ── MESH → NLM Linked Data API (D-number → label) ────────────────────────
    if "mesh" in src:
        try:
            r = requests.get(f"https://id.nlm.nih.gov/mesh/{eid}.json", timeout=8)
            if r.ok:
                data = r.json()
                label = (data.get("prefLabel")
                         or (data.get("label") or {}).get("@value", ""))
                if label:
                    return entity_id, label
        except Exception:
            pass
        # Fallback: NCBI eSearch by descriptor UI
        try:
            sr = requests.get(f"{NCBI_BASE}/esearch.fcgi",
                params={"db": "mesh", "term": f"{eid}[MH]",
                        "retmax": 1, "retmode": "json"}, timeout=6)
            uid_list = sr.json().get("esearchresult", {}).get("idlist", [])
            if uid_list:
                smr = requests.get(f"{NCBI_BASE}/esummary.fcgi",
                    params={"db": "mesh", "id": uid_list[0], "retmode": "json"},
                    timeout=6)
                info = smr.json().get("result", {}).get(uid_list[0], {})
                name = info.get("ds_meshui_name") or info.get("ds_name", "")
                if name:
                    return entity_id, name
        except Exception:
            pass
        return entity_id, eid

    # ── OMIM / other: use the ID as label ────────────────────────────────────
    return entity_id, eid


def _resolve_disease_names(entity_ids: List[str],
                            max_workers: int = 8) -> Dict[str, str]:
    """Parallel name resolution for a list of disease entity IDs."""
    results: Dict[str, str] = {}
    total = len(entity_ids)
    done  = 0

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(_resolve_one_disease, eid): eid
                   for eid in entity_ids}
        for fut in as_completed(futures):
            eid, name = fut.result()
            results[eid] = name
            done += 1
            if done % 50 == 0 or done == total:
                print(f"  Resolved {done}/{total} disease names…")

    return results


# ── Phase 3: neuro keyword filter ─────────────────────────────────────────────

def _is_neuro(name: str) -> bool:
    n = name.lower()
    return any(kw in n for kw in NEURO_KEYWORDS)


# ── Compound + gene name resolution (unchanged logic) ─────────────────────────

_DBID_RE = re.compile(r"^DB\d{5,}$")


def _is_drugbank(entity: str) -> bool:
    return "drugbank" in entity.lower() or _DBID_RE.match(_eid(entity)) is not None


def _resolve_compounds(cmp_ids: List[str]) -> Dict[str, dict]:
    result: Dict[str, dict] = {}
    db_ids = [c for c in cmp_ids if _is_drugbank(c)]
    other  = [c for c in cmp_ids if not _is_drugbank(c)]

    print(f"  Resolving {len(db_ids)} DrugBank compound names via ChEMBL…")
    for i, eid in enumerate(db_ids):
        db_id = _eid(eid)
        try:
            r = requests.get(f"{CHEMBL_BASE}/molecule.json",
                params={"cross_references__xref_src": "DrugBank",
                        "cross_references__xref_id": db_id,
                        "limit": 1, "format": "json"}, timeout=8)
            if r.ok:
                mols = r.json().get("molecules", [])
                if mols:
                    m = mols[0]
                    result[eid] = {
                        "name":       (m.get("pref_name") or db_id).upper(),
                        "chembl_id":  m.get("molecule_chembl_id", ""),
                        "drugbank_id": db_id,
                        "max_phase":  m.get("max_phase") or 0,
                    }
                    continue
        except Exception:
            pass
        result[eid] = {"name": db_id, "chembl_id": "",
                        "drugbank_id": db_id, "max_phase": 0}
        if (i + 1) % 50 == 0:
            print(f"    {i+1}/{len(db_ids)}…")
        time.sleep(0.12)

    for eid in other:
        result[eid] = {"name": _eid(eid), "chembl_id": "",
                        "drugbank_id": "", "max_phase": 0}
    return result


def _resolve_genes(gene_ids: List[str]) -> Dict[str, dict]:
    result: Dict[str, dict] = {}
    entrez = [g for g in gene_ids if "entrez" in g.lower()]

    print(f"  Resolving {len(entrez)} gene symbols via NCBI…")
    for i in range(0, len(entrez), 200):
        batch   = entrez[i:i + 200]
        id_nums = [_eid(g) for g in batch]
        try:
            r = requests.get(f"{NCBI_BASE}/esummary.fcgi",
                params={"db": "gene", "id": ",".join(id_nums), "retmode": "json"},
                timeout=20)
            if r.ok:
                data = r.json().get("result", {})
                for eid in batch:
                    enum = _eid(eid)
                    info = data.get(enum, {})
                    result[eid] = {
                        "symbol":    info.get("name", "") or enum,
                        "gene_name": info.get("description", ""),
                        "entrez_id": enum,
                    }
        except Exception:
            for eid in batch:
                result[eid] = {"symbol": _eid(eid), "gene_name": "", "entrez_id": _eid(eid)}
        time.sleep(0.35)

    for eid in gene_ids:
        if eid not in result:
            result[eid] = {"symbol": _eid(eid), "gene_name": "", "entrez_id": ""}
    return result


# ── Public API ─────────────────────────────────────────────────────────────────

def load_cache() -> dict:
    try:
        if CACHE_FILE.exists():
            return json.loads(CACHE_FILE.read_text(encoding="utf-8"))
    except Exception:
        pass
    return {}


def cache_available() -> bool:
    return CACHE_FILE.exists() and CACHE_FILE.stat().st_size > 10_000


def build_cache(max_compounds: int = 400, local_file: str = None) -> dict:
    """
    Dynamically build DRKG neuro subset cache.
    1. Parse all 5.8M triplets → discover every disease entity
    2. Resolve disease names via OLS (DOID) / NCBI (MESH)
    3. Filter to neurological diseases by keyword
    4. Keep compound + gene edges for those diseases
    5. Resolve compound (DrugBank→ChEMBL) and gene (Entrez→symbol) names
    """
    print("=" * 64)
    print("DRKG  ·  Dynamic Neurological Disease Subset Builder")
    print("=" * 64)

    # ── Load raw triplets ──────────────────────────────────────────────────────
    if local_file:
        p = Path(local_file).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(f"Not found: {p}")
        lines = _load_tsv_from_file(p)
    else:
        lines = _download_tsv()

    # ── Phase 1: discover entities ────────────────────────────────────────────
    print(f"\nPhase 1 — Scanning {len(lines):,} triplets…")
    dis_count, treat_all, dg_all, target_all, all_cmp_ids = \
        _discover_disease_entities(lines)

    print(f"  Unique disease entities found : {len(dis_count):,}")
    print(f"  Unique compound entities      : {len(all_cmp_ids):,}")
    print(f"  Compound-disease edges (raw)  : {len(treat_all):,}")
    print(f"  Disease-gene edges (raw)      : {len(dg_all):,}")

    # ── Phase 2: resolve disease names (parallel) ─────────────────────────────
    # Only resolve diseases that actually appear in edges (performance)
    # Focus on top 1000 by edge count — they're the most connected/important
    top_dis_ids = [eid for eid, _ in dis_count.most_common(1000)]

    print(f"\nPhase 2 — Resolving names for {len(top_dis_ids)} disease entities…")
    dis_name_map = _resolve_disease_names(top_dis_ids, max_workers=8)

    # ── Phase 3: filter to neuro diseases ─────────────────────────────────────
    neuro_ids: Set[str] = {
        eid for eid, name in dis_name_map.items() if _is_neuro(name)
    }

    print(f"\nPhase 3 -- Neurological diseases identified: {len(neuro_ids)}")
    for eid in sorted(neuro_ids):
        print(f"  {eid}  ->  {dis_name_map[eid]}")

    if not neuro_ids:
        print("\n  WARNING: 0 neuro diseases matched. Check keyword list.")
        print("  Sample entity IDs found:")
        for eid, cnt in dis_count.most_common(10):
            print(f"    {eid}  (count={cnt})")
        return {}

    # ── Phase 4: filter edges to neuro diseases ───────────────────────────────
    treat_edges = [e for e in treat_all if e["t"] in neuro_ids]
    dg_edges    = [e for e in dg_all    if e["h"] in neuro_ids]

    # Count compounds per disease for ranking
    cmp_degree: Counter = Counter(e["h"] for e in treat_edges)

    # Cap to top compounds
    top_cmps    = {eid for eid, _ in cmp_degree.most_common(max_compounds)}
    treat_edges = [e for e in treat_edges if e["h"] in top_cmps]
    target_edges = [e for e in target_all if e["h"] in top_cmps]

    # Final gene set
    gene_ids: Set[str] = set()
    for e in dg_edges:    gene_ids.add(e["t"])
    for e in target_edges: gene_ids.add(e["t"])

    print(f"\n  Treat edges   : {len(treat_edges):,}")
    print(f"  D-G edges     : {len(dg_edges):,}")
    print(f"  Target edges  : {len(target_edges):,}")
    print(f"  Compounds     : {len(top_cmps)}")
    print(f"  Genes         : {len(gene_ids)}")

    # ── Phase 5: resolve compound + gene names ────────────────────────────────
    print("\nPhase 5 — Resolving compound and gene names…")
    compounds = _resolve_compounds(list(top_cmps))
    genes     = _resolve_genes(list(gene_ids))

    # ── Build disease info dict ───────────────────────────────────────────────
    diseases = {
        eid: {"name": dis_name_map.get(eid, _eid(eid)), "entity_id": eid,
              "id_type": _esrc(eid), "id_value": _eid(eid)}
        for eid in neuro_ids
    }

    data = {
        "diseases":        diseases,
        "compounds":       compounds,
        "genes":           genes,
        "edges": {
            "treat":    treat_edges,
            "dis_gene": dg_edges,
            "target":   target_edges,
        },
        "disease_name_index": {
            info["name"].lower(): eid for eid, info in diseases.items()
        },
        "_ts": time.time(),
        "_counts": {
            "diseases":    len(diseases),
            "compounds":   len(compounds),
            "genes":       len(genes),
            "treat_edges": len(treat_edges),
            "dg_edges":    len(dg_edges),
            "tgt_edges":   len(target_edges),
        },
    }

    CACHE_FILE.parent.mkdir(exist_ok=True)
    CACHE_FILE.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")

    print("\nCache saved to", CACHE_FILE)
    for k, v in data["_counts"].items():
        print(f"   {k:<20}: {v:,}")

    return data


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    local = sys.argv[1] if len(sys.argv) > 1 else None
    build_cache(local_file=local)
