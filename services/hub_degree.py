"""
Hub-degree normalization  —  DWPC-style down-weighting of network hubs.
═══════════════════════════════════════════════════════════════════════════════
A network HUB (EGFR deg~506, TP53~766) sits in dozens of pathways and interacts
with hundreds of proteins, so it trivially maxes out the pathway and PPI cohesion
scores for EVERY disease it touches — non-discriminating inflation, not evidence.
The Rephetio / Himmelstein (2017) fix is degree-weighting: a gene's contribution
along a path is damped by its connectivity, `degree^(-w)` with w=0.4.

Here we apply it to the LIVE pathway/PPI scores via a per-gene damping factor,
normalized so a typical-degree gene (~REF) keeps full weight and a hub is damped:

    damping(gene) = (REF / max(REF, global_PPI_degree))^0.4

Global degree is the high-confidence STRING PPI partner count (cached). FAIL-SOFT:
unknown degree → damping 1.0 (no change), so a data gap never invents a penalty.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

import requests

logger = logging.getLogger(__name__)

STRING_PARTNERS = "https://string-db.org/api/json/interaction_partners"
_REF_DEGREE = 50.0        # a typical (non-hub) gene degree — full weight at/below this
_DAMP = 0.4               # Himmelstein 2017 degree-damping exponent
_deg_cache: Dict[str, Optional[int]] = {}


def gene_degree(gene: str) -> Optional[int]:
    """Global high-confidence PPI degree of a gene (STRING partner count). Cached;
    None on failure (→ no damping)."""
    g = (gene or "").strip().upper()
    if not g:
        return None
    if g in _deg_cache:
        return _deg_cache[g]
    deg = None
    try:
        r = requests.get(STRING_PARTNERS, params={
            "identifiers": g, "species": 9606, "required_score": 700,
            "limit": 2000, "caller_identity": "neurorepurpose_platform"}, timeout=12)
        if r.ok:
            deg = len(r.json())
    except Exception as e:
        logger.debug(f"STRING degree failed for {gene}: {e}")
    _deg_cache[g] = deg
    return deg


def damping(gene: str) -> float:
    """Degree-damping factor in (0,1]: 1.0 for a typical/specific gene, lower for a
    hub. Fail-soft → 1.0 when degree is unknown."""
    deg = gene_degree(gene)
    if deg is None or deg <= _REF_DEGREE:
        return 1.0
    return round((_REF_DEGREE / deg) ** _DAMP, 3)


def damping_map(genes: List[str]) -> Dict[str, float]:
    """{GENE_UPPER: damping} for a set of genes (typically a drug's few targets)."""
    return {g.upper(): damping(g) for g in (genes or [])}
