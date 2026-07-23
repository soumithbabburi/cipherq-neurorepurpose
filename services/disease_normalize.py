"""
Disease normalization + subtype expansion (platform-wide).
═══════════════════════════════════════════════════════════════════════════════
The rule, applied everywhere a disease is the input:

    Normalize the term to a controlled ontology, then split a broad PARENT label
    into its specific SUBTYPES before scoring or reporting candidates.

Broad umbrella terms ("vasculitis", "cancer", "arthritis") hide clinically distinct
subtypes (vasculitis → EGPA, GPA, MPA, Takayasu, GCA, Behçet). Scoring/ranking a
single broad label is imprecise and misleading. This module detects a broad term
and expands it to subtype-level indications, keeping the parent only as a category.

Backbone: MeSH provides the parent/child hierarchy we already hold (mesh_diseases),
MONDO/EFO the canonical id (via disease_ontology / disease_id). Subtypes are ranked
by the Repurposing Value Score so the most pursued indications surface first.

Fail-soft: no DB / unresolved term → treated as a specific (non-broad) indication,
so callers behave exactly as before.
"""
from __future__ import annotations

import logging
import re
from functools import lru_cache
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# A term is treated as broad when it has at least this many ontology children.
_MIN_CHILDREN_BROAD = 4
# How many subtypes to expand into (keeps candidate counts + latency bounded).
_MAX_SUBTYPES = 8

# Curated umbrella terms that are clinically too broad to score as one indication —
# a safety net for well-known categories even if the MeSH child count is modest.
# Matched as whole terms (not substrings) so "breast cancer" is NOT caught by "cancer".
_BROAD_TERMS = {
    "vasculitis", "cancer", "neoplasm", "neoplasms", "solid tumor", "arthritis",
    "dementia", "leukemia", "lymphoma", "anemia", "neuropathy", "myopathy",
    "muscular dystrophy", "inflammatory bowel disease", "colitis", "dermatitis",
    "nephritis", "hepatitis", "encephalitis", "neurodegenerative disease",
    "autoimmune disease", "interstitial lung disease", "cardiomyopathy",
    "epilepsy", "movement disorder", "connective tissue disease",
}


def _norm_text(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "").strip().lower())


def _children(mesh_id: str) -> List[Dict]:
    """Direct ontology children (subtypes) of a MeSH disease id."""
    try:
        from services.disease_resolver import _q
    except Exception:
        return []
    if not mesh_id or str(mesh_id).startswith(("LOCAL:", "MONDO", "EFO")):
        return []
    rows = _q(
        "SELECT mesh_id, heading, tree_numbers FROM mesh_diseases "
        "WHERE parent_ids && ARRAY[%s]::text[] ORDER BY heading",
        (mesh_id,),
    )
    # drop self-referential / duplicate headings
    seen, out = set(), []
    for r in rows:
        h = (r.get("heading") or "").strip()
        if h and h.lower() not in seen:
            seen.add(h.lower())
            out.append(r)
    return out


def _value_of(name: str) -> float:
    """Repurposing Value Score (0–1) for ranking subtypes; 0 when unavailable."""
    try:
        from services.disease_value import value_for
        v = value_for(name) or {}
        return float(v.get("value_score") or 0.0)
    except Exception:
        return 0.0


@lru_cache(maxsize=2048)
def normalize(name: str) -> Dict:
    """Resolve a disease name to its canonical record + broadness signal.

    Returns {input, label, mesh_id, is_broad, n_children}. `is_broad` is True when the
    term is an umbrella category (many ontology children, or a curated umbrella term).
    """
    q = _norm_text(name)
    out = {"input": name, "label": name.strip(), "mesh_id": "", "mondo_id": "",
           "n_children": 0, "is_broad": False, "tree_numbers": ()}
    if not q:
        return out
    try:
        from services.disease_resolver import resolve_disease as _mesh_resolve
        recs = _mesh_resolve(name) or []
    except Exception as e:
        logger.debug("normalize: mesh resolve failed: %s", e)
        recs = []
    rec = recs[0] if recs else {}
    if rec.get("unresolved"):
        rec = {}
    label = (rec.get("heading") or name).strip()
    mesh_id = rec.get("mesh_id", "") or ""
    # MONDO backbone (complete + orphan-free) — primary hierarchy when built.
    mondo_id, mondo_children = "", 0
    try:
        from services import mondo_hierarchy as _mh
        if _mh.available():
            mondo_id = _mh.resolve(label, mesh_id) or _mh.resolve(name, mesh_id) or ""
            if mondo_id:
                mondo_children = _mh.child_count(mondo_id)
    except Exception as e:
        logger.debug("normalize: MONDO resolve failed: %s", e)
    n_children = len(_children(mesh_id))     # MeSH child count (fallback signal)
    # Broad = an umbrella category (many MONDO or MeSH children) or a curated umbrella
    # term. A specific term like "breast cancer" (which has children) is still expandable
    # (→ TNBC/HER2+/…); a true leaf ("EGPA") has ~0 children and stays as-is.
    is_broad = (mondo_children >= _MIN_CHILDREN_BROAD) or (n_children >= _MIN_CHILDREN_BROAD) \
        or (_norm_text(label) in _BROAD_TERMS) or (q in _BROAD_TERMS)
    out.update(label=label, mesh_id=mesh_id, mondo_id=mondo_id,
               n_children=max(n_children, mondo_children), is_broad=bool(is_broad),
               tree_numbers=tuple(rec.get("tree_numbers") or ()))
    return out


@lru_cache(maxsize=2048)
def subtypes(name: str, limit: int = _MAX_SUBTYPES) -> tuple:
    """Clinically-distinct subtypes of a (broad) disease, ranked by Repurposing Value.
    Returns a tuple of {label, mesh_id, value_score} (tuple so it is lru-cacheable)."""
    info = normalize(name)
    # PRIMARY: MONDO hierarchy (complete + orphan-free — reaches EGPA under vasculitis,
    # which MeSH's orphan row cannot). Falls through to the MeSH tree only if MONDO has
    # no id / no subtypes for this term.
    if info.get("mondo_id"):
        try:
            from services import mondo_hierarchy as _mh
            ms = _mh.subtypes(info["mondo_id"], limit)
            if ms:
                return tuple({"label": s["label"], "mesh_id": s.get("mesh_id", ""),
                              "value_score": s["value_score"]} for s in ms)
        except Exception as e:
            logger.debug("subtypes: MONDO path failed: %s", e)

    if not info["mesh_id"] or not info.get("tree_numbers"):
        return tuple()
    # FALLBACK — descend via MeSH TREE NUMBERS (more complete than parent_ids, which has
    # orphan rows): every descendant of vasculitis (C14.907.940) shares that prefix, so
    # one query pulls GCA/GPA/MPA/Takayasu/Behçet. We then keep the LEAVES — a node whose
    # tree number is not a prefix of another descendant's — i.e. the clinically-distinct
    # subtypes, dropping intermediate umbrella nodes (Systemic Vasculitis, ANCA-associated).
    try:
        from services.disease_resolver import _q
    except Exception:
        return tuple()
    descendants: Dict[str, Dict] = {}
    all_tns: set = set()
    for tn in info["tree_numbers"]:
        rows = _q("SELECT mesh_id, heading, tree_numbers FROM mesh_diseases "
                  "WHERE EXISTS (SELECT 1 FROM unnest(tree_numbers) t WHERE t LIKE %s)",
                  (tn + ".%",))
        for r in rows:
            h = (r.get("heading") or "").strip()
            if not h:
                continue
            descendants.setdefault(h.lower(), r)
            for t in (r.get("tree_numbers") or []):
                all_tns.add(t)
    # a descendant is a LEAF if no other descendant tree number extends it
    def _is_leaf(r) -> bool:
        for t in (r.get("tree_numbers") or []):
            if any(o != t and o.startswith(t + ".") for o in all_tns):
                return False
        return True
    leaves = [r for r in descendants.values() if _is_leaf(r)]
    scored = [{"label": (r.get("heading") or "").strip(), "mesh_id": r.get("mesh_id", ""),
               "value_score": round(_value_of((r.get("heading") or "").strip()), 3)}
              for r in leaves if (r.get("heading") or "").strip()]
    # rank by repurposing value (most pursued indications first), then name for stability
    scored.sort(key=lambda s: (-s["value_score"], s["label"]))
    return tuple(scored[:limit])


def expand(name: str, limit: int = _MAX_SUBTYPES) -> Dict:
    """The platform-wide entry point. Returns the set of INDICATIONS to actually score:

        {is_broad, parent_label, parent_mesh_id,
         indications: [{label, mesh_id, is_subtype}]}

    For a specific term the indications list is just [the term itself]; for a broad
    umbrella it is the ranked subtypes (parent kept only as `parent_label`). Callers
    iterate `indications` and score/report each, grouping under `parent_label`.
    """
    info = normalize(name)
    if not info["is_broad"]:
        return {"is_broad": False, "parent_label": info["label"],
                "parent_mesh_id": info["mesh_id"],
                "indications": [{"label": info["label"], "mesh_id": info["mesh_id"],
                                 "is_subtype": False}]}
    subs = subtypes(name, limit)
    if not subs:                       # broad but no resolvable children — degrade safely
        return {"is_broad": False, "parent_label": info["label"],
                "parent_mesh_id": info["mesh_id"],
                "indications": [{"label": info["label"], "mesh_id": info["mesh_id"],
                                 "is_subtype": False}]}
    return {
        "is_broad": True,
        "parent_label": info["label"],
        "parent_mesh_id": info["mesh_id"],
        "indications": [{"label": s["label"], "mesh_id": s["mesh_id"],
                         "is_subtype": True, "value_score": s["value_score"]} for s in subs],
    }


if __name__ == "__main__":
    for d in ["vasculitis", "cancer", "breast cancer", "rheumatoid arthritis",
              "eosinophilic granulomatosis with polyangiitis", "leukemia"]:
        e = expand(d)
        tag = "BROAD" if e["is_broad"] else "specific"
        print(f"\n{d}  [{tag}]  parent={e['parent_label']}")
        for ind in e["indications"][:8]:
            print(f"   - {ind['label']}" + (f"  (value {ind.get('value_score')})" if ind.get('is_subtype') else ""))
