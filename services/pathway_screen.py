"""
Pathway-first repurposing engine
═══════════════════════════════════════════════════════════════════════════════
Lets a user enter through a BIOLOGICAL PATHWAY instead of a disease, and still
returns genuine repurposing candidates — i.e. existing drug → NEW indication
pairs, with the pathway as the mechanistic rationale.

A pathway by itself is not a repurposing result; repurposing always needs a
drug↔indication pair. So this engine treats the pathway as a *bridge*:

  Mode B  (pathway-first discovery, default)
  ─────────────────────────────────────────
    pathway ─▶ member genes ─▶ drugs that modulate those genes (with direction)
            └▶ member genes ─▶ diseases the pathway drives (Open Targets)
    The repurposing claim = (drug, inferred disease) the drug is NOT already
    approved for, scored with the SAME canonical pair score the rest of the
    platform uses.

  Mode A  (disease-anchored refinement)
  ─────────────────────────────────────
    Same, but the caller supplies the disease. The pathway then acts purely as
    a mechanistic filter/ranker on a normal repurposing screen for that disease.

Everything here is additive and reuses existing platform machinery:
  • data/protein_pathways.json + data/pathways.json  → offline pathway↔gene index
  • services.repurposing_engine._get_chembl_pool      → local chembl_33 (drugs/action)
  • services.signature_engine                          → signed modulation direction
  • services.reverse_repurposing                       → gene→disease + canonical pair score
  • services.regulatory_verdict / services.pos_model   → per-candidate enrichment
"""
from __future__ import annotations

import json
import logging
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Bounded fan-out for the network-bound enrichment loops (Open Targets / ChEMBL).
# These calls are I/O-bound, so threads collapse tens of seconds of serial waiting
# into a few seconds without overwhelming the upstream APIs.
_MAX_WORKERS = 8
# The (drug × disease) pair scorer is the heaviest loop; give it a wider pool since
# each task is almost entirely network wait, then cap how many inferred diseases we
# score per drug (they're pre-ranked, so the best claim is near the top).
_MAX_PAIR_WORKERS = 14
_MAX_ENRICH_DISEASES = 5

logger = logging.getLogger(__name__)

_DATA = Path(__file__).parent.parent / "data"
PATHWAYS_FILE  = Path(__file__).parent.parent / "pathways.json"
PROTPATH_FILE  = Path(__file__).parent.parent / "protein_pathways.json"
CACHE_FILE     = _DATA / "pathway_cache.json"
CACHE_TTL      = 21600  # 6h, matching the other engines


# ── Screen cache ───────────────────────────────────────────────────────────────

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


# ── Pathway ↔ gene index (offline, built once) ─────────────────────────────────
# pathways.json     : "Pathway::PC7_1063" -> {"name": "...", "source": "Reactome", ...}
# protein_pathways  : "GENE" -> ["Pathway::PC7_1063", ...]
# We invert protein_pathways into pathway_id -> [genes], and keep a name index.

_pw_names: Dict[str, dict] = {}          # pathway_id -> meta
_pw_genes: Dict[str, List[str]] = {}     # pathway_id -> [gene symbols]
_name_index: List[Tuple[str, str]] = []  # (lowercased name, pathway_id)
_norm_index: List[Tuple[frozenset, str]] = []  # (normalized token set, pathway_id)
_loaded = False

# ── Fuzzy name normalization ────────────────────────────────────────────────
# Curated pathway names spell the same concept many ways (NF-kappaB / NF-kB / NFkB,
# TGF-beta / TGFB, signaling / signalling). We reduce a name to a canonical token
# SET so a user's spelling resolves to the pathway regardless of these variants.
# Greek letters are collapsed to a single latin letter on BOTH the query and the
# index so "NF-kB" and "NF-kappaB" both become {nf, kb}.
_GREEK_UNICODE = {"α": "a", "β": "b", "γ": "g", "δ": "d",
                  "κ": "k", "λ": "l", "σ": "s", "ω": "w"}
_GREEK_WORDS = [("kappa", "k"), ("alpha", "a"), ("beta", "b"), ("gamma", "g"),
                ("delta", "d"), ("lambda", "l"), ("sigma", "s"), ("omega", "w")]
# Filler tokens carry no discriminating power for matching (spelling variants of
# "signaling" included), so they're dropped from the token set on both sides.
_PW_FILLER = {"signaling", "signalling", "signal", "signals", "pathway", "pathways",
              "the", "of", "via", "by", "in", "and", "a", "an", "to"}


def _pw_norm_tokens(name: str) -> frozenset:
    s = (name or "").lower()
    for u, l in _GREEK_UNICODE.items():
        s = s.replace(u, l)
    for w, l in _GREEK_WORDS:
        s = s.replace(w, l)
    toks = re.split(r"[^a-z0-9]+", s)
    return frozenset(t for t in toks if t and t not in _PW_FILLER)


def _ensure_index() -> bool:
    global _loaded
    if _loaded:
        return bool(_pw_genes)
    try:
        names = json.loads(PATHWAYS_FILE.read_text(encoding="utf-8"))
        prot  = json.loads(PROTPATH_FILE.read_text(encoding="utf-8"))
    except Exception as e:
        logger.warning(f"pathway index load failed: {e}")
        _loaded = True
        return False

    for pid, meta in names.items():
        nm = (meta.get("name") or "").strip()
        _pw_names[pid] = {"name": nm, "source": meta.get("source", ""),
                          "category": meta.get("category", "")}
        if nm:
            _name_index.append((nm.lower(), pid))
            _norm_index.append((_pw_norm_tokens(nm), pid))

    for gene, pids in prot.items():
        g = (gene or "").strip().upper()
        if not g:
            continue
        for pid in pids:
            _pw_genes.setdefault(pid, []).append(g)

    # Dedupe + sort gene lists for stable output
    for pid in list(_pw_genes):
        _pw_genes[pid] = sorted(set(_pw_genes[pid]))

    _loaded = True
    logger.info(f"pathway index: {len(_pw_names)} pathways, "
                f"{sum(len(v) for v in _pw_genes.values())} gene-memberships")
    return bool(_pw_genes)


def _short_id(pid: str) -> str:
    """'Pathway::PC7_1063' -> 'PC7_1063' for display."""
    return pid.split("::", 1)[-1] if "::" in pid else pid


def _fuzzy_pathway(query: str) -> Optional[str]:
    """Best pathway id for a query that failed exact/substring matching.

    Reduces the query to a canonical token set and finds pathways whose token set
    CONTAINS all of them (spelling/greek/synonym tolerant). Ties broken toward the
    tightest match (fewest extra tokens) and then the most member genes. Falls back
    to strong partial overlap (Jaccard) when no full-containment candidate exists."""
    qt = _pw_norm_tokens(query)
    if not qt:
        return None

    contained: List[Tuple[int, int, str]] = []   # (extra_tokens, -genes, pid)
    partial:   List[Tuple[float, int, str]] = []  # (-jaccard, -genes, pid)
    for toks, pid in _norm_index:
        if not toks:
            continue
        ngenes = len(_pw_genes.get(pid, []))
        if qt <= toks:                       # every query token present
            contained.append((len(toks - qt), -ngenes, pid))
        else:
            inter = len(qt & toks)
            if inter:
                jac = inter / len(qt | toks)
                # require a meaningful overlap so we don't match on one shared token
                if jac >= 0.5 or (inter >= 2 and inter >= len(qt) - 1):
                    partial.append((-jac, -ngenes, pid))

    if contained:
        contained.sort()
        return contained[0][2]
    if partial:
        partial.sort()
        return partial[0][2]
    return None


def suggest_pathways(prefix: str, limit: int = 12) -> List[Dict]:
    """Autocomplete: pathways whose name contains the query, ranked by gene count."""
    if not _ensure_index():
        return []
    q = (prefix or "").strip().lower()
    if not q:
        return []
    hits = []
    for nm, pid in _name_index:
        if q in nm:
            hits.append({"id": _short_id(pid), "full_id": pid,
                         "name": _pw_names[pid]["name"],
                         "n_genes": len(_pw_genes.get(pid, []))})
    hits.sort(key=lambda x: (-(q == x["name"].lower()), -x["n_genes"]))

    # Fuzzy fallback (greek / spelling / synonym) when substring finds nothing —
    # keeps autocomplete useful for queries like "NF-kB signaling".
    if not hits:
        qt = _pw_norm_tokens(prefix)
        if qt:
            fuzzy = []
            for toks, pid in _norm_index:
                if toks and qt <= toks:
                    fuzzy.append((len(toks - qt), -len(_pw_genes.get(pid, [])), pid))
            fuzzy.sort()
            for _, _, pid in fuzzy[:limit]:
                hits.append({"id": _short_id(pid), "full_id": pid,
                             "name": _pw_names[pid]["name"],
                             "n_genes": len(_pw_genes.get(pid, []))})
    return hits[:limit]


def resolve_pathway(query: str) -> Dict:
    """Resolve a pathway NAME or ID to {id, full_id, name, genes:[...], n_genes}.

    Accepts a full id ('Pathway::PC7_1063'), a short id ('PC7_1063'), or a name
    (exact, then substring). Returns {} when nothing matches with member genes."""
    if not _ensure_index():
        return {}
    q = (query or "").strip()
    if not q:
        return {}

    pid = None
    # 1. direct id forms
    if q in _pw_names:
        pid = q
    elif f"Pathway::{q}" in _pw_names:
        pid = f"Pathway::{q}"
    else:
        ql = q.lower()
        # 2. exact name
        for nm, candidate in _name_index:
            if nm == ql:
                pid = candidate
                break
        # 3. substring name — pick the matching pathway with the most genes
        if pid is None:
            subs = [c for nm, c in _name_index if ql in nm]
            if subs:
                pid = max(subs, key=lambda c: len(_pw_genes.get(c, [])))

        # 4. fuzzy token match — spelling / greek / synonym tolerance.
        #    e.g. "NF-kB signaling" -> {nf, kb} ⊆ "Canonical NF-kappaB pathway".
        if pid is None:
            pid = _fuzzy_pathway(q)

    if not pid:
        return {}
    genes = _pw_genes.get(pid, [])
    if not genes:
        return {}
    return {"id": _short_id(pid), "full_id": pid,
            "name": _pw_names[pid]["name"], "source": _pw_names[pid].get("source", ""),
            "genes": genes, "n_genes": len(genes)}


# ── Pathway genes → drugs that modulate them (local chembl_33, REST fallback) ──

def _drugs_for_genes(genes: List[str]) -> List[Dict]:
    """Every drug with a curated mechanism on any of `genes`, with action direction.

    Returns rows: {chembl_id, name, max_phase, gene, action_type, mechanism}.
    One row per (drug, gene, action). Local chembl_33 first; ChEMBL REST fallback."""
    upper = sorted({g.strip().upper() for g in genes if g and g.strip()})
    if not upper:
        return []

    # ── Local chembl_33 (fast, offline) ──────────────────────────────────────
    try:
        from services.repurposing_engine import _get_chembl_pool
        pool = _get_chembl_pool()
    except Exception:
        pool = None

    if pool is not None:
        conn = None
        try:
            conn = pool.getconn()
            with conn.cursor() as cur:
                cur.execute(
                    """
                    SELECT md.chembl_id, md.pref_name, md.max_phase,
                           UPPER(csyn.component_synonym) AS gene,
                           dm.action_type, dm.mechanism_of_action
                    FROM component_synonyms csyn
                    JOIN target_components tc ON tc.component_id = csyn.component_id
                    JOIN drug_mechanism dm    ON dm.tid = tc.tid
                    JOIN molecule_dictionary md ON md.molregno = dm.molregno
                    WHERE csyn.syn_type = 'GENE_SYMBOL'
                      AND UPPER(csyn.component_synonym) = ANY(%s)
                      AND md.pref_name IS NOT NULL
                    """,
                    (upper,),
                )
                rows = [
                    {"chembl_id": r[0], "name": r[1] or r[0],
                     "max_phase": float(r[2]) if r[2] is not None else 0.0,
                     "gene": r[3], "action_type": r[4] or "",
                     "mechanism": r[5] or ""}
                    for r in cur.fetchall()
                ]
            if rows:
                return rows
        except Exception as e:
            logger.debug(f"local drugs-for-genes failed: {e}")
        finally:
            if conn is not None:
                pool.putconn(conn)

    # ── REST fallback: gene → target → mechanisms → molecules ─────────────────
    return _drugs_for_genes_rest(upper)


def _drugs_for_genes_rest(upper: List[str]) -> List[Dict]:
    """ChEMBL REST fallback for _drugs_for_genes (bounded to a handful of genes)."""
    from services import http_client
    base = "https://www.ebi.ac.uk/chembl/api/data"
    out: List[Dict] = []
    for gene in upper[:20]:
        try:
            tr = http_client.get(f"{base}/target.json",
                                 params={"target_synonym__iexact": gene, "limit": 5,
                                         "format": "json"}, timeout=10)
            tids = [t.get("target_chembl_id") for t in
                    (tr.json().get("targets", []) if tr and tr.ok else [])
                    if t.get("target_chembl_id")]
            if not tids:
                continue
            mr = http_client.get(f"{base}/mechanism.json",
                                 params={"target_chembl_id__in": ",".join(tids[:5]),
                                         "limit": 200, "format": "json"}, timeout=12)
            mechs = mr.json().get("mechanisms", []) if mr and mr.ok else []
            mids = sorted({m.get("molecule_chembl_id") for m in mechs
                           if m.get("molecule_chembl_id")})
            if not mids:
                continue
            from services.repurposing_engine import _chembl_molecule_details
            details = _chembl_molecule_details(mids)
            for m in mechs:
                mid = m.get("molecule_chembl_id")
                d = details.get(mid, {})
                if not mid or not d.get("name"):
                    continue
                out.append({"chembl_id": mid, "name": d.get("name", mid),
                            "max_phase": float(d.get("max_phase") or 0),
                            "gene": gene, "action_type": m.get("action_type", "") or "",
                            "mechanism": m.get("mechanism_of_action", "") or ""})
        except Exception as e:
            logger.debug(f"REST drugs-for-gene {gene} failed: {e}")
    return out


# ── Direction handling ─────────────────────────────────────────────────────────

def _signed(action: str) -> float:
    from services.signature_engine import _direction
    return _direction(action)


def _net_direction(signed_vals: List[float]) -> str:
    s = sum(signed_vals)
    if s < -1e-9:
        return "suppresses"
    if s > 1e-9:
        return "activates"
    return "mixed"


def _direction_match(net: str, desired: str) -> float:
    """How well a drug's net effect matches what the caller wants on the pathway."""
    if desired == "either":
        return 1.0 if net in ("suppresses", "activates") else 0.4
    want = "suppresses" if desired == "suppress" else "activates"
    if net == want:
        return 1.0
    if net == "mixed":
        return 0.4
    return 0.0  # opposes the desired direction


# ── Pathway → inferred diseases (mode B) ───────────────────────────────────────

# Umbrella terms that are too generic to be an actionable repurposing indication —
# Open Targets associates them strongly with almost any oncogenic/inflammatory
# pathway, but "repurpose drug X for cancer" is not a claim a client can file.
_GENERIC_DISEASE = {
    "cancer", "neoplasm", "neoplasms", "tumor", "tumour", "carcinoma",
    "malignant neoplasm", "disease", "disorder", "syndrome",
    "cell proliferation disorder", "hematologic cancer", "hematopoietic neoplasm",
    "immune system disease", "nervous system disease", "genetic disease",
}


def _is_generic_disease(name: str) -> bool:
    n = (name or "").strip().lower()
    return n in _GENERIC_DISEASE or len(n) < 4


def _inferred_diseases(genes: List[str], max_genes: int = 15,
                       top: int = 8) -> List[Dict]:
    """Diseases the pathway drives: aggregate Open Targets target→disease
    associations across the pathway's member genes. Cached per gene via the
    reverse engine's ontology cache."""
    try:
        from services.reverse_repurposing import _gene_to_ensembl, _diseases_for_target
    except Exception as e:
        logger.debug(f"inferred-diseases import failed: {e}")
        return []

    def _rows_for_gene(gene: str):
        try:
            ens = _gene_to_ensembl(gene)
            if not ens:
                return gene, []
            return gene, _diseases_for_target(ens, size=100)
        except Exception as e:
            logger.debug(f"inferred-diseases gene {gene} failed: {e}")
            return gene, []

    target_genes = genes[:max_genes]
    agg: Dict[str, Dict] = {}
    with ThreadPoolExecutor(max_workers=min(_MAX_WORKERS, len(target_genes) or 1)) as ex:
        for gene, rows in ex.map(_rows_for_gene, target_genes):
            for row in rows:
                name = row.get("name", "")
                if not name or _is_generic_disease(name):
                    continue
                key = name.lower()
                e = agg.setdefault(key, {
                    "name": name, "efo": row.get("efo_id", ""),
                    "areas": row.get("therapeutic_areas", []) or [],
                    "score_sum": 0.0, "genes": []})
                e["score_sum"] += float(row.get("score", 0) or 0)
                if gene not in e["genes"]:
                    e["genes"].append(gene)

    # Attach the Repurposing Value Score so a pathway's inferred indications are ranked
    # by "would a pharma company pursue this disease?" (burden × unmet-need × market),
    # not only by how many pathway genes drive it — the same filter the reverse screen
    # uses, so trivial inferred indications don't top the pathway results either.
    try:
        from services.disease_value import value_for
    except Exception:
        value_for = None
    for e in agg.values():
        dv = value_for(e["name"], e.get("efo", "")) if value_for else None
        e["disease_value"] = dv
        e["_vw"] = (dv or {}).get("value_score", 0.5) if dv else 0.5
    ranked = sorted(agg.values(),
                    key=lambda x: (x["_vw"] * (0.4 + 0.15 * len(x["genes"])), x["score_sum"]),
                    reverse=True)
    for r in ranked:
        r["score_sum"] = round(r["score_sum"], 3)
        r["n_pathway_genes"] = len(r["genes"])
        r.pop("_vw", None)
    return ranked[:top]


# ── Main engine ────────────────────────────────────────────────────────────────

def screen_pathway(query: str, direction: str = "either",
                   disease: Optional[str] = None,
                   max_drugs: int = 40, enrich_top: int = 12) -> Dict:
    """Pathway-first repurposing screen.

    direction : 'suppress' | 'activate' | 'either'  (desired effect on the pathway)
    disease   : optional — supplying it switches to Mode A (disease-anchored).
    """
    direction = (direction or "either").lower().strip()
    if direction not in ("suppress", "activate", "either"):
        direction = "either"
    mode = "A" if disease else "B"

    cache_key = f"pwscreen:{query.lower().strip()}|{direction}|{(disease or '').lower().strip()}"
    cache = _load_cache()
    if cache_key in cache:
        return cache[cache_key]

    pathway = resolve_pathway(query)
    if not pathway:
        return {"error": f"Pathway not recognised: '{query}'. Try a name like "
                         f"'mTOR signaling' or a pathway id.", "candidates": []}

    genes = pathway["genes"]

    # ── Candidate drugs that modulate the pathway, aggregated per drug ─────────
    rows = _drugs_for_genes(genes)
    by_drug: Dict[str, Dict] = {}
    for r in rows:
        d = by_drug.setdefault(r["chembl_id"], {
            "chembl_id": r["chembl_id"], "name": r["name"],
            "max_phase": r["max_phase"], "nodes_hit": [], "_seen": set()})
        sig = _signed(r["action_type"])
        node = {"gene": r["gene"], "action": r["action_type"] or "modulator",
                "mechanism": r["mechanism"],
                "direction": ("down" if sig < 0 else "up" if sig > 0 else "neutral")}
        key = (r["gene"], r["action_type"])
        if key not in d["_seen"]:
            d["_seen"].add(key)
            d["nodes_hit"].append(node)

    if not by_drug:
        result = {"pathway": pathway, "direction": direction, "mode": mode,
                  "candidates": [], "inferred_diseases": [],
                  "error": "No drugs with curated mechanisms on this pathway's genes."}
        return result

    pw_gene_set = set(genes)

    # ── Score modulation (mechanistic coverage + direction match) ─────────────
    cands: List[Dict] = []
    for cid, d in by_drug.items():
        d.pop("_seen", None)
        nodes = d["nodes_hit"]
        signed_vals = [_signed(n["action"]) for n in nodes]
        net = _net_direction(signed_vals)
        n_hit = len({n["gene"] for n in nodes})
        coverage = n_hit / max(1, len(pw_gene_set))
        # coverage saturates (hitting 3 nodes is already strong evidence of pathway action)
        cov_score = 1.0 - pow(2.71828, -n_hit / 2.0)
        dmatch = _direction_match(net, direction)
        phase_w = min(1.0, float(d["max_phase"] or 0) / 4.0)
        modulation = round(0.55 * dmatch + 0.30 * cov_score + 0.15 * phase_w, 4)
        d.update({
            "n_nodes_hit": n_hit,
            "coverage": round(coverage, 3),
            "net_direction": net,
            "direction_match": round(dmatch, 2),
            "modulation_score": modulation,
        })
        cands.append(d)

    # Drop drugs that clearly oppose the desired direction (when a direction is set)
    if direction != "either":
        cands = [c for c in cands if c["direction_match"] > 0.0]

    cands.sort(key=lambda x: (x["modulation_score"], x["n_nodes_hit"],
                              x["max_phase"]), reverse=True)
    cands = cands[:max_drugs]

    # ── Resolve the indication space (the part that makes this REPURPOSING) ────
    inferred = [] if mode == "A" else _inferred_diseases(genes)

    # In disease-anchored mode, understand a manifestation-qualified anchor (e.g.
    # "systemic sclerosis with pulmonary manifestations") as its PRIMARY disorder,
    # so the whole screen (and the score) is anchored on the treatable disease and
    # the UI can label it a partial fit rather than silently mis-scoring the phrase.
    anchor_resolution: Dict = {}
    if mode == "A" and disease:
        try:
            from services.disease_ontology import resolve_disease as _resolve_disease
            di = _resolve_disease(disease) or {}
            anchor_resolution = {
                "queried_as": disease,
                "resolved_name": di.get("disease_name", ""),
                "partial_fit": bool(di.get("partial_fit")),
                "primary_disorder": di.get("primary_disorder", ""),
                "manifestation": di.get("manifestation", ""),
                "classification": di.get("classification", {}),
            }
        except Exception as e:
            logger.debug(f"anchor resolution failed: {e}")

    # ── Attach the repurposing claim (drug → NEW indication) ──────────────────
    _attach_indications(cands[:enrich_top], mode, disease, inferred)

    result = {
        "pathway": pathway,
        "direction": direction,
        "mode": mode,
        "anchor_disease": disease or "",
        "anchor_resolution": anchor_resolution,
        "n_candidates": len(cands),
        "inferred_diseases": inferred,
        "candidates": cands,
        "method": ("Pathway-first repurposing: pathway→genes→drugs (direction-aware) "
                   "with indications inferred from Open Targets gene→disease "
                   "associations; pairs scored by the canonical repurposing score."),
        "_ts": time.time(),
    }
    cache[cache_key] = result
    _save_cache(cache)
    return result


def _attach_indications(top: List[Dict], mode: str, disease: Optional[str],
                        inferred: List[Dict]):
    """For each top drug, attach its repurposing claim + enrichment.

    Mode A : score against the supplied disease (skip if already approved there).
    Mode B : pick the best inferred disease the drug is NOT already approved for.
    """
    try:
        from services.reverse_repurposing import canonical_pair_score, resolve_drug
    except Exception as e:
        logger.debug(f"indication attach import failed: {e}")
        return

    try:
        from services.regulatory_verdict import assess as reg_assess, _matches
    except Exception:
        reg_assess = None
        def _matches(a, b): return (a or "").lower() in (b or "").lower()
    try:
        from services.pos_model import predict_progression
    except Exception:
        predict_progression = None

    if not top:
        return

    # The candidate indications are already ranked upstream; scoring every drug
    # against all of them is the dominant cost. Cap to the strongest few — the
    # best novel claim virtually always sits among the top-ranked diseases.
    targets_disease = ([disease] if mode == "A"
                       else [d["name"] for d in inferred])[:_MAX_ENRICH_DISEASES]
    targets_disease = [d for d in targets_disease if d]

    # ── Phase 1: resolve each drug once (parallel) ────────────────────────────
    def _resolve(c):
        try:
            return c["chembl_id"], (resolve_drug(c["chembl_id"]) or {})
        except Exception:
            return c["chembl_id"], {}

    infos: Dict[str, Dict] = {}
    with ThreadPoolExecutor(max_workers=min(_MAX_WORKERS, len(top))) as ex:
        for cid, info in ex.map(_resolve, top):
            infos[cid] = info

    # ── Phase 2: pre-warm each unique disease's context once ──────────────────
    # canonical_pair_score resolves the disease (Open Targets + STRING) internally;
    # warming it once per disease stops N concurrent pairs from re-fetching it.
    if targets_disease:
        try:
            from services.disease_ontology import resolve_disease
            with ThreadPoolExecutor(max_workers=min(_MAX_WORKERS, len(targets_disease))) as ex:
                list(ex.map(lambda d: _safe(resolve_disease, d), targets_disease))
        except Exception:
            pass

    # ── Phase 3: score every (drug, disease) pair in one flat parallel pool ────
    def _score(pair):
        c, dis = pair
        cid = c["chembl_id"]
        drug_genes = sorted({n["gene"] for n in c["nodes_hit"]})
        known_list = [k.get("name", "") for k in
                      (infos.get(cid, {}).get("known_indications", []) or []) if k.get("name")]
        try:
            ps = canonical_pair_score(
                chembl_id=cid, disease=dis, drug_genes=drug_genes,
                max_phase=int(float(c.get("max_phase") or 0)),
                indications="; ".join(known_list), drug_name=c.get("name", ""),
                # The pathway screen has ALREADY established that this drug
                # direction-matches a pathway driving `dis`; pass that modulation as
                # the mechanistic prior so the pair score credits the linkage that
                # surfaced it instead of scoring it as zero-mechanism.
                mechanistic_prior=c.get("modulation_score"))
            # Genuine repurposing requires a NEW use: flag indications the drug is
            # already developed/approved for (mechanistic confirmation, not a lead).
            is_known = any(_matches(dis, k) for k in known_list)
            return cid, {"disease": dis,
                         "repurposing_score": ps.get("score", 0.0),
                         "scores": ps.get("scores", {}),
                         "is_novel": not is_known}
        except Exception as e:
            logger.debug(f"pair score {cid}/{dis} failed: {e}")
            return cid, None

    pairs = [(c, dis) for c in top for dis in targets_disease]
    by_cid: Dict[str, List[Dict]] = {c["chembl_id"]: [] for c in top}
    if pairs:
        with ThreadPoolExecutor(max_workers=min(_MAX_PAIR_WORKERS, len(pairs))) as ex:
            for cid, row in ex.map(_score, pairs):
                if row:
                    by_cid[cid].append(row)

    # ── Phase 4: aggregate + local enrichment (serial, cheap) ─────────────────
    for c in top:
        cid = c["chembl_id"]
        # Rank novel indications first, then by score — the repurposing claim is the
        # top NOVEL indication; known ones are kept but marked as confirmatory.
        scored_inds = sorted(by_cid.get(cid, []),
                             key=lambda x: (x["is_novel"], x["repurposing_score"]),
                             reverse=True)
        c["candidate_indications"] = scored_inds[:3]
        novel = [s for s in scored_inds if s["is_novel"]]
        if novel:
            c["best_indication"] = novel[0]["disease"]
            c["repurposing_score"] = novel[0]["repurposing_score"]
            c["best_is_novel"] = True
        elif scored_inds:
            c["best_indication"] = scored_inds[0]["disease"]
            c["repurposing_score"] = scored_inds[0]["repurposing_score"]
            c["best_is_novel"] = False

        known = infos.get(cid, {}).get("known_indications", []) or []
        if reg_assess:
            try:
                c["regulatory"] = reg_assess(
                    c.get("name", ""), c.get("max_phase", 0), known,
                    c.get("best_indication", ""))
            except Exception:
                pass
        if predict_progression:
            try:
                c["pos"] = predict_progression(
                    current_phase=float(c.get("max_phase") or 0),
                    evidence_score=c.get("repurposing_score"),
                    is_repurposing=True)
            except Exception:
                pass
        # Same platform quality filters — potency/lead-viability + the value of the
        # best inferred indication — so the pathway surface is consistent with the rest.
        if c.get("best_indication"):
            try:
                from services.quality_overlay import overlay
                c.update(overlay(c.get("name", ""), cid, c["best_indication"],
                                 sorted({n["gene"] for n in c.get("nodes_hit", [])})))
            except Exception:
                pass


def _safe(fn, *args):
    try:
        return fn(*args)
    except Exception:
        return None


if __name__ == "__main__":
    import sys
    q = sys.argv[1] if len(sys.argv) > 1 else "mTOR signaling"
    d = sys.argv[2] if len(sys.argv) > 2 else "suppress"
    res = screen_pathway(q, direction=d)
    if res.get("error"):
        print("ERROR:", res["error"])
    else:
        pw = res["pathway"]
        print(f"Pathway: {pw['name']} ({pw['id']}) — {pw['n_genes']} genes, mode {res['mode']}")
        print(f"Inferred diseases: {[d['name'] for d in res['inferred_diseases']]}")
        print(f"Top candidates ({res['n_candidates']}):")
        for c in res["candidates"][:10]:
            print(f"  {c['name']:<28} mod={c['modulation_score']:.2f} "
                  f"{c['net_direction']:<11} nodes={c['n_nodes_hit']} "
                  f"→ {c.get('best_indication','?')} ({c.get('repurposing_score','-')})")
