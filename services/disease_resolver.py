"""
Disease Resolver
Accepts any free-text disease query and resolves it to a canonical MeSH ID
plus all related IDs (parents, children, siblings) using the mesh_diseases table.
Falls back to LIKE-based matching if MeSH lookup fails.
"""

import logging
import os
import re
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

import psycopg2
import psycopg2.extras
import psycopg2.pool

logger = logging.getLogger(__name__)

from config import db_params  # centralized DB config (no hardcoded credentials)

_pool: Optional[psycopg2.pool.ThreadedConnectionPool] = None
_pool_failed = False


def _get_pool():
    global _pool, _pool_failed
    if _pool is None:
        if _pool_failed:
            return None
        try:
            _pool = psycopg2.pool.ThreadedConnectionPool(1, 8, **db_params())   # headroom for the parallel reverse-screen scoring loop
        except Exception as e:
            _pool_failed = True
            logger.warning(f"DB pool failed: {e}")
    return _pool


def _q(sql: str, params=None) -> List[Dict]:
    pool = _get_pool()
    if pool is None:
        return []
    conn = pool.getconn()
    try:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(sql, params)
            return [dict(r) for r in cur.fetchall()]
    except Exception as e:
        logger.error(f"Query failed: {e}")
        return []
    finally:
        pool.putconn(conn)


# Common abbreviations and alternate names → canonical form
_ALIASES = {
    "ad": "alzheimer disease",
    "alzheimers": "alzheimer disease",
    "alzheimer's": "alzheimer disease",
    "pd": "parkinson disease",
    "parkinsons": "parkinson disease",
    "parkinson's": "parkinson disease",
    "ms": "multiple sclerosis",
    "als": "amyotrophic lateral sclerosis",
    "lou gehrig": "amyotrophic lateral sclerosis",
    "mnd": "motor neuron disease",
    "hd": "huntington disease",
    "huntingtons": "huntington disease",
    "lbd": "lewy body disease",
    "ftd": "frontotemporal dementia",
    "mci": "mild cognitive impairment",
    "bpd": "bipolar disorder",
    "ocd": "obsessive-compulsive disorder",
    "ptsd": "stress disorders, post-traumatic",
    "adhd": "attention deficit disorder with hyperactivity",
    "asd": "autism spectrum disorder",
    "tbi": "brain injuries, traumatic",
    "cte": "chronic traumatic encephalopathy",
    # Non-neuro research areas (oncology / cardiometabolic / autoimmune / respiratory / rare)
    "t2d": "diabetes mellitus, type 2",
    "t2dm": "diabetes mellitus, type 2",
    "ra": "arthritis, rheumatoid",
    "sle": "lupus erythematosus, systemic",
    "lupus": "lupus erythematosus, systemic",
    "nsclc": "carcinoma, non-small-cell lung",
    "aml": "leukemia, myeloid, acute",
    "crc": "colorectal neoplasms",
    "hcc": "carcinoma, hepatocellular",
    "gbm": "glioblastoma",
    "mm": "multiple myeloma",
    "copd": "pulmonary disease, chronic obstructive",
    "ipf": "idiopathic pulmonary fibrosis",
    "cf": "cystic fibrosis",
    "sma": "muscular atrophy, spinal",
    "dmd": "muscular dystrophy, duchenne",
    "nash": "non-alcoholic fatty liver disease",
    "nafld": "non-alcoholic fatty liver disease",
    "ckd": "renal insufficiency, chronic",
    "pah": "hypertension, pulmonary",
    "afib": "atrial fibrillation",
    "cad": "coronary artery disease",
    "nmo": "neuromyelitis optica",
}


def _normalize(query: str) -> str:
    q = query.strip().lower()
    q = re.sub(r"['’]s?\b", "", q)  # remove possessive
    q = re.sub(r"\s+", " ", q).strip()
    return _ALIASES.get(q, q)


# ── Primary-disorder / systemic-manifestation parsing ───────────────────────
# Clinicians phrase an indication as a PRIMARY (often systemic) disorder qualified
# by the organ/system it manifests in — e.g. "systemic sclerosis with pulmonary
# manifestations", "amyloidosis with cardiac involvement", "lupus-related nephritis".
# The disease databases key on the primary disorder, so a literal lookup of the whole
# phrase misses. We extract the primary disorder and treat it as a PARTIAL FIT: the
# rationale still points at the primary disorder, but the caller is told it was a
# partial (manifestation-qualified) match, not an exact indication.
_MANIFEST_WORDS = (r"manifestation|involvement|feature|symptom|complication|sign|"
                   r"presentation|sequela|phenotype|damage|disease|disorder")
_PRIMARY_PATTERNS = [
    # "<primary> with <manifestation> manifestations/involvement/..."
    re.compile(rf"^(?P<primary>.+?)\s+with\s+(?P<manifestation>.+?)\s+(?:{_MANIFEST_WORDS})s?\b",
               re.I),
    # "<manifestation> manifestations/involvement of/in <primary>"
    re.compile(rf"^(?P<manifestation>.+?)\s+(?:{_MANIFEST_WORDS})s?\s+(?:of|in|due to)\s+(?P<primary>.+)$",
               re.I),
    # "<primary> with <organ/system> involvement"  (involvement without a middle noun)
    re.compile(r"^(?P<primary>.+?)\s+with\s+(?P<manifestation>.+?)\s+involvement$", re.I),
    # "<primary>-related|-associated|-driven <manifestation>" — the systemic driver
    # (e.g. lupus) is the primary disorder; the trailing noun is the manifestation.
    re.compile(r"^(?P<primary>.+?)[- ](?:related|associated|driven)\s+(?P<manifestation>.+)$",
               re.I),
    # "<primary> secondary to <manifestation>"  (primary is the disorder we can treat)
    re.compile(r"^(?P<primary>.+?)\s+secondary to\s+(?P<manifestation>.+)$", re.I),
]


def parse_primary_disorder(query: str) -> Optional[Tuple[str, str]]:
    """Split a manifestation-qualified disease phrase into (primary_disorder,
    manifestation). Returns None when the query is a single, unqualified concept.

    'systemic sclerosis with pulmonary manifestations' -> ('systemic sclerosis', 'pulmonary')
    'cardiac involvement in amyloidosis'               -> ('amyloidosis', 'cardiac')
    'lupus-related nephritis'                          -> ('lupus', 'nephritis')  (via alias -> SLE)
    'alzheimer disease'                                -> None
    """
    q = (query or "").strip()
    if not q or len(q) < 6:
        return None
    for pat in _PRIMARY_PATTERNS:
        m = pat.search(q)
        if not m:
            continue
        primary = (m.groupdict().get("primary") or "").strip(" ,-.")
        manifestation = (m.groupdict().get("manifestation") or "").strip(" ,-.")
        # The primary must be a substantive disorder name, and the split must have
        # actually removed something (otherwise it's the same concept).
        if primary and len(primary) >= 4 and primary.lower() != q.lower():
            return primary, manifestation
    return None


# ── Primary organ-specific vs multi-systemic classification ─────────────────
# The platform must distinguish:
#   • Primary Organ-Specific Pathology — core etiology and tissue destruction
#     originate within a SINGLE target organ system (e.g. vitiligo → skin,
#     idiopathic pulmonary fibrosis → lung, glaucoma → eye, ulcerative colitis → gut).
#   • Multi-Systemic / Syndromic Disorder with Secondary Organ Manifestations —
#     local findings are downstream of a BROADER systemic disease (e.g. SLE,
#     sarcoidosis, systemic sclerosis, amyloidosis).
# The signal is data-driven, not a hand-kept disease list: MeSH tree numbers place
# each disease in one or more of the C-category organ systems. One localizable
# system ⇒ organ-specific; several systems, or an inherently-systemic branch
# (connective-tissue, immune, blood, metabolic), ⇒ multi-systemic. Manifestation
# phrasing and a few name markers refine borderline calls.

# MeSH C-category → human organ-system label
_MESH_SYSTEM = {
    "C01": "infectious", "C04": "neoplastic", "C05": "musculoskeletal",
    "C06": "digestive", "C07": "oral / dental", "C08": "respiratory",
    "C09": "ear-nose-throat", "C10": "nervous system", "C11": "eye",
    "C12": "urinary / male reproductive", "C13": "female reproductive",
    "C14": "cardiovascular", "C15": "blood & lymphatic", "C16": "congenital",
    "C17": "skin / connective tissue", "C18": "metabolic", "C19": "endocrine",
    "C20": "immune", "C23": "general pathology", "C25": "chemically-induced",
    "C26": "injury",
}
# Localizable systems: a disease confined here is organ-specific.
_LOCALIZABLE = {"C06", "C07", "C08", "C09", "C10", "C11", "C12", "C13", "C14"}
# Branches that are systemic by nature even when they appear alone.
_INHERENT_SYSTEMIC = {"C15": "blood & lymphatic", "C18": "metabolic", "C20": "immune"}
# Dominant-organ label for the organ-specific verdict.
_TARGET_ORGAN = {
    "C06": "gastrointestinal tract", "C07": "mouth / dentition", "C08": "lungs",
    "C09": "ear-nose-throat", "C10": "nervous system", "C11": "eye",
    "C12": "genitourinary tract", "C13": "female reproductive tract",
    "C14": "heart & vasculature", "C19": "endocrine gland", "skin": "skin",
}
_SYSTEMIC_NAME = ("systemic", "syndrome", "disseminated", "generalized", "generalised",
                  "diffuse", "multiorgan", "multi-organ", "multisystem", "vasculitis",
                  "polyangiitis", "connective tissue", "sarcoid", "amyloid", "lupus")


def _systems_from_trees(tree_numbers: List[str]):
    """Split MeSH tree numbers into (localizable organ systems, systemic-branch flags)."""
    systems, systemic = set(), set()
    for tn in tree_numbers or []:
        top = tn.split(".")[0]
        if top == "C17":                       # skin (C17.800) vs connective tissue (C17.300)
            if tn.startswith("C17.300"):
                systemic.add("connective tissue")
            elif tn.startswith("C17.800"):
                systems.add("skin")
            else:
                systems.add("skin / connective tissue")
        elif top in _INHERENT_SYSTEMIC:
            systemic.add(_INHERENT_SYSTEMIC[top])
        elif top in _LOCALIZABLE:
            systems.add(_MESH_SYSTEM.get(top, top))
        elif top == "C19":                     # single endocrine gland = localizable
            systems.add("endocrine")
        elif top in _MESH_SYSTEM:
            systemic.add(_MESH_SYSTEM[top])
    return systems, systemic


def classify_disease(query: str, rows: Optional[List[Dict]] = None) -> Dict:
    """Classify a disease as Primary Organ-Specific vs Multi-Systemic / Syndromic.

    Returns {category, confidence, target_organ, affected_systems, systemic_driver,
    manifestation, rationale, signals}. Category is one of:
        'primary_organ_specific' | 'multi_systemic_syndromic' | 'indeterminate'.
    """
    q = (query or "").strip()
    if rows is None:
        rows = resolve_disease(q)
    rec = rows[0] if rows else {}
    heading = rec.get("heading", q)
    trees = rec.get("tree_numbers") or []
    resolved = bool(trees) and not rec.get("unresolved")

    signals: List[str] = []

    # Signal 1 — manifestation phrasing (a systemic disorder qualified by an organ)
    pf = parse_primary_disorder(q)
    manifestation = ""
    systemic_driver = ""
    if pf:
        systemic_driver, manifestation = pf
        signals.append(f"phrased as a disorder with '{manifestation}' manifestations")

    # Signal 2 — name markers
    name_l = f"{q} {heading}".lower()
    name_systemic = any(k in name_l for k in _SYSTEMIC_NAME)
    if name_systemic:
        signals.append("name carries a systemic / syndromic marker")

    # Signal 3 — MeSH organ-system breadth
    systems, systemic_branches = _systems_from_trees(trees)
    n_local = len(systems)
    if trees:
        touched = sorted(systems | systemic_branches)
        signals.append(f"MeSH places it in: {', '.join(touched) or 'unclassified'}")

    multi = bool(pf) or name_systemic or bool(systemic_branches) or n_local >= 2
    affected = sorted(systems | systemic_branches)

    if multi:
        category = "multi_systemic_syndromic"
        # best-effort dominant organ of the secondary manifestation
        target_organ = (next(iter(systems)) if len(systems) == 1 else "")
        confidence = "high" if (resolved and (bool(systemic_branches) or n_local >= 2)) \
            else ("medium" if (resolved or pf or name_systemic) else "low")
        rationale = ("Multi-systemic / syndromic: local findings are downstream of a broader "
                     "systemic disease" +
                     (f" ({systemic_driver})" if systemic_driver else "") +
                     (f"; secondary organ system(s): {', '.join(affected)}" if affected else "") + ".")
    elif n_local == 1:
        category = "primary_organ_specific"
        organ = next(iter(systems))
        # map back to a target-organ noun where we can
        target_organ = organ
        for code, sysname in _MESH_SYSTEM.items():
            if sysname == organ and code in _TARGET_ORGAN:
                target_organ = _TARGET_ORGAN[code]
                break
        if organ == "skin":
            target_organ = "skin"
        confidence = "high" if resolved else "low"
        rationale = (f"Primary organ-specific: core etiology and tissue destruction originate "
                     f"within the {target_organ} ({organ} system); no systemic driver detected.")
    else:
        category = "indeterminate"
        target_organ = ""
        confidence = "low"
        rationale = ("Insufficient ontology signal to classify organ-specific vs systemic; "
                     "treated as indeterminate.")

    return {
        "category": category,
        "label": {"primary_organ_specific": "Primary Organ-Specific",
                  "multi_systemic_syndromic": "Multi-Systemic / Syndromic",
                  "indeterminate": "Indeterminate"}[category],
        "confidence": confidence,
        "heading": heading,
        "target_organ": target_organ,
        "affected_systems": affected,
        "systemic_driver": systemic_driver,
        "manifestation": manifestation,
        "rationale": rationale,
        "signals": signals,
    }


def resolve_disease(query: str) -> List[Dict]:
    """
    Resolve a free-text query to one or more MeSH disease records.
    Returns list of dicts with: mesh_id, heading, tree_numbers, entry_terms, parent_ids
    """
    norm = _normalize(query)

    # 1. Exact heading match
    rows = _q(
        "SELECT * FROM mesh_diseases WHERE LOWER(heading) = %s",
        (norm,),
    )
    if rows:
        return rows

    # 2. Exact entry term match
    rows = _q(
        "SELECT * FROM mesh_diseases WHERE %s = ANY(LOWER(entry_terms::text)::text[])",
        (norm,),
    )
    if not rows:
        rows = _q(
            "SELECT * FROM mesh_diseases WHERE EXISTS ("
            "  SELECT 1 FROM unnest(entry_terms) t WHERE LOWER(t) = %s"
            ")",
            (norm,),
        )
    if rows:
        return rows

    # 3. Partial heading match
    rows = _q(
        "SELECT * FROM mesh_diseases WHERE LOWER(heading) LIKE %s ORDER BY heading LIMIT 10",
        (f"%{norm}%",),
    )
    if rows:
        return rows

    # 4. Partial entry term match
    rows = _q(
        "SELECT * FROM mesh_diseases WHERE EXISTS ("
        "  SELECT 1 FROM unnest(entry_terms) t WHERE LOWER(t) LIKE %s"
        ") ORDER BY heading LIMIT 10",
        (f"%{norm}%",),
    )
    if rows:
        return rows

    # 5. Partial fit — a primary/systemic disorder qualified by manifestations.
    #    Resolve the primary disorder and flag the result as a partial (not exact) fit
    #    so downstream scoring/UI can label it honestly.
    pf = parse_primary_disorder(query)
    if pf:
        primary, manifestation = pf
        sub = resolve_disease(primary)
        if sub and not sub[0].get("unresolved"):
            out = []
            for r in sub:
                r = dict(r)
                r["partial_fit"] = True
                r["primary_disorder"] = primary
                r["manifestation"] = manifestation
                r["query"] = query.strip()
                out.append(r)
            return out

    # 6. Word-by-word fallback — match first significant word
    words = [w for w in norm.split() if len(w) > 3]
    if words:
        pattern = f"%{words[0]}%"
        rows = _q(
            "SELECT * FROM mesh_diseases WHERE LOWER(heading) LIKE %s ORDER BY heading LIMIT 10",
            (pattern,),
        )
        if rows:
            return rows

    return [{
        "mesh_id": norm,
        "heading": query.strip() or norm,
        "tree_numbers": [],
        "entry_terms": [],
        "parent_ids": [],
        "unresolved": True,
    }]


def expand_mesh_ids(mesh_ids: List[str], include_children: bool = True) -> List[str]:
    """
    Given a list of MeSH IDs, return all related IDs:
    - The IDs themselves
    - Their parents (broader disease category)
    - Their children (more specific subtypes) — if include_children=True
    """
    if not mesh_ids:
        return []

    if all(str(mid).startswith("LOCAL:") for mid in mesh_ids):
        return mesh_ids

    all_ids = set(mesh_ids)

    # Add parents (one level up)
    for mid in mesh_ids:
        rows = _q("SELECT parent_ids FROM mesh_diseases WHERE mesh_id = %s", (mid,))
        for r in rows:
            for pid in (r.get("parent_ids") or []):
                all_ids.add(pid)

    if include_children:
        # Add children: any disease whose parent_ids contains one of our IDs
        placeholders = ",".join(["%s"] * len(mesh_ids))
        rows = _q(
            f"SELECT mesh_id FROM mesh_diseases "
            f"WHERE parent_ids && ARRAY[{placeholders}]::text[]",
            mesh_ids,
        )
        for r in rows:
            all_ids.add(r["mesh_id"])

    return list(all_ids)


def get_tree_siblings(mesh_id: str) -> List[str]:
    """Return MeSH IDs that share the same immediate parent (siblings in the tree)."""
    rows = _q("SELECT parent_ids FROM mesh_diseases WHERE mesh_id = %s", (mesh_id,))
    if not rows or not rows[0].get("parent_ids"):
        return []
    parent_ids = rows[0]["parent_ids"]
    if not parent_ids:
        return []
    placeholders = ",".join(["%s"] * len(parent_ids))
    sibling_rows = _q(
        f"SELECT mesh_id FROM mesh_diseases "
        f"WHERE parent_ids && ARRAY[{placeholders}]::text[] AND mesh_id != %s",
        parent_ids + [mesh_id],
    )
    return [r["mesh_id"] for r in sibling_rows]


def get_disease_label(mesh_id: str) -> str:
    """Return the canonical heading for a MeSH ID."""
    if str(mesh_id).startswith("LOCAL:"):
        return str(mesh_id).replace("LOCAL:", "").replace("-", " ").title()

    rows = _q("SELECT heading FROM mesh_diseases WHERE mesh_id = %s", (mesh_id,))
    return rows[0]["heading"] if rows else mesh_id


def suggest_diseases(partial: str, limit: int = 10) -> List[str]:
    """Return disease name suggestions for autocomplete."""
    if not partial or len(partial) < 2:
        return []
    norm = partial.strip().lower()
    rows = _q(
        "SELECT DISTINCT heading FROM mesh_diseases "
        "WHERE LOWER(heading) LIKE %s ORDER BY heading LIMIT %s",
        (f"%{norm}%", limit),
    )
    suggestions = [r["heading"] for r in rows]
    if suggestions:
        return suggestions
    fallback = [
        "Alzheimer Disease",
        "Parkinson Disease",
        "Multiple Sclerosis",
        "Epilepsy",
        "Depression",
    ]
    return [s for s in fallback if norm in s.lower()][:limit]


def mesh_available() -> bool:
    """Check whether the mesh_diseases table has data."""
    rows = _q("SELECT COUNT(*) AS n FROM mesh_diseases")
    return bool(rows) and rows[0]["n"] > 0
