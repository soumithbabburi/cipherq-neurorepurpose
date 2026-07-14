"""
Target-driven pathogenic proliferation  —  a mechanism the genetic-overlap scorer is
structurally blind to.
═══════════════════════════════════════════════════════════════════════════════
Some repurposing works not because the drug's target is a GENETIC-ASSOCIATION gene
for the disease, but because the target DRIVES a pathogenic cellular process the
disease is built on — hyperproliferation. Palbociclib → rheumatoid arthritis is the
canonical miss: CDK4/6–cyclin-D1 drives synovial fibroblast (FLS) hyperplasia, so a
CDK4/6 INHIBITOR is a real lead — yet CDK4/6 isn't an RA "gene", so target overlap = 0
and the pair gets phantom-gated to ~0.01.

This module recognises that mechanism, DATA-DRIVEN and DIRECTION-AWARE, generalising to
any cell-cycle inhibitor into any hyperproliferative disease (RA, fibrosis, restenosis,
psoriasis, systemic sclerosis, glomerulonephritis, BPH …):

  proliferation_match fires ONLY when   (a) the drug's target is a POSITIVE proliferation
  driver (an accelerator you INHIBIT to suppress growth — not a tumour suppressor),
  AND (b) the drug's action on it is INHIBITORY, AND (c) the disease is characterised by
  pathogenic hyperproliferation.

Signals (all cached, fail-soft):
  • gene is a cell-cycle gene   — local Reactome cell-cycle pathway membership
  • gene DIRECTION              — textbook activator/suppressor sets (GO's own
                                  proliferation-regulation terms are too noisy — GO even
                                  labels CDK6 a "negative" regulator), corroborated by GO
  • disease hyperproliferation  — text (name/def/synonyms/OT description) + a PubMed
                                  literature ratio (RA→8500 vs Alzheimer→14)
"""
from __future__ import annotations

import json
import logging
import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

_ROOT = Path(__file__).parent.parent


# ── Gene side: is the target a POSITIVE proliferation driver you'd inhibit? ──────

# Textbook core cell-cycle ACCELERATORS — inhibiting them suppresses proliferation.
# (Direction is stable biology; used as the GO-noise-proof backstop, not the detector.)
_ACTIVATORS = {
    "CDK1", "CDK2", "CDK3", "CDK4", "CDK6", "CDK7", "CDK9",
    "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2", "CCNA1", "CCNA2", "CCNB1", "CCNB2",
    "E2F1", "E2F2", "E2F3", "MYC", "MYCN", "MYCL", "FOXM1",
    "AURKA", "AURKB", "PLK1", "PLK4", "BUB1", "TTK", "WEE1",
    "CDC25A", "CDC25B", "CDC25C", "CDC20", "CDC6", "CDC7", "SKP2", "MDM2",
    "CDK4/6", "MKI67", "TOP2A", "PCNA", "RRM2",
}
# Tumour SUPPRESSORS / negative regulators — inhibiting them would PROMOTE proliferation.
_SUPPRESSORS = {
    "RB1", "RBL1", "RBL2", "TP53", "TP63", "TP73", "PTEN",
    "CDKN1A", "CDKN1B", "CDKN1C", "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D",
    "CHEK1", "CHEK2", "ATM", "ATR", "APC", "NF1", "NF2", "VHL", "WT1",
}
# Actions (ChEMBL action_type) that suppress a target's activity.
_INHIBITORY = ("inhibitor", "antagonist", "blocker", "negative modulator",
               "inverse agonist", "degrader", "disruptor")


@lru_cache(maxsize=1)
def _cell_cycle_genes() -> frozenset:
    """Genes in a core cell-cycle Reactome pathway (dynamic membership, local data)."""
    try:
        pw = json.loads((_ROOT / "pathways.json").read_text(encoding="utf-8"))
        prot = json.loads((_ROOT / "protein_pathways.json").read_text(encoding="utf-8"))
    except Exception as e:
        logger.debug(f"cell-cycle gene load failed: {e}")
        return frozenset()
    cc_ids = {k for k, v in pw.items()
              if re.search(r"cell cycle|cell division|g1/s|mitotic|dna replication",
                           (v.get("name") or "").lower())}
    out = {g.upper() for g, pids in prot.items() if any(p in cc_ids for p in pids)}
    return frozenset(out)


@lru_cache(maxsize=2048)
def _go_direction(gene: str) -> Optional[str]:
    """'positive' / 'negative' / None from MyGene GO proliferation terms (corroboration
    only — noisy, so never overrides the textbook sets)."""
    try:
        import requests
        r = requests.get("https://mygene.info/v3/query",
                         params={"q": gene, "species": "human", "fields": "go", "size": 1},
                         timeout=12)
        hits = r.json().get("hits", []) if r.ok else []
        go = hits[0].get("go", {}) if hits else {}
        ids = set()
        for aspect in ("BP", "MF", "CC"):
            v = go.get(aspect)
            for t in (v if isinstance(v, list) else [v] if v else []):
                if t.get("id"):
                    ids.add(t["id"])
        if "GO:0008284" in ids or "GO:0045787" in ids:   # positive reg proliferation / cell cycle
            return "positive"
        if "GO:0008285" in ids or "GO:0045786" in ids or "GO:0007050" in ids:  # negative / arrest
            return "negative"
    except Exception as e:
        logger.debug(f"GO direction {gene} failed: {e}")
    return None


def proliferation_role(gene: str) -> Dict:
    """{is_driver, direction, confidence, evidence} for a gene. A 'positive driver' is an
    accelerator whose inhibition suppresses proliferation (CDK4/6, cyclins, MYC…)."""
    g = (gene or "").strip().upper()
    out = {"gene": g, "is_driver": False, "direction": "unknown",
           "confidence": 0.0, "evidence": []}
    if not g:
        return out
    in_reactome = g in _cell_cycle_genes()
    if g in _ACTIVATORS:
        out.update(is_driver=True, direction="positive", confidence=0.95,
                   evidence=["core cell-cycle activator (textbook)"])
        if in_reactome:
            out["evidence"].append("Reactome cell-cycle pathway")
        return out
    if g in _SUPPRESSORS:
        out.update(is_driver=False, direction="negative", confidence=0.9,
                   evidence=["tumour suppressor / negative regulator (textbook)"])
        return out
    # Not textbook — allow a DYNAMIC positive call only with agreeing Reactome + GO,
    # so we never miscredit an unknown gene that happens to be a suppressor.
    if in_reactome:
        go = _go_direction(g)
        if go == "positive":
            out.update(is_driver=True, direction="positive", confidence=0.6,
                       evidence=["Reactome cell-cycle pathway", "GO positive-proliferation"])
        elif go == "negative":
            out.update(direction="negative", confidence=0.6,
                       evidence=["Reactome cell-cycle pathway", "GO negative/arrest"])
        else:
            out.update(direction="unknown", confidence=0.3,
                       evidence=["Reactome cell-cycle pathway (direction unresolved)"])
    return out


# ── Disease side: is the disease characterised by pathogenic hyperproliferation? ─

_PROLIF_KW = (
    "prolifera", "hyperprolifera", "hyperplasia", "hyperplastic", "fibros",
    "fibrotic", "pannus", "synovial", "sclero", "stenosis", "restenosis",
    "neointima", "hypertroph", "keloid", "granulation tissue", "dysplasia",
    "angiogen", "remodel", "overgrowth", "smooth muscle", "keratinocyte",
    "mesangial", "myofibroblast", "granuloma",
)


def _text_signal(*texts: str) -> List[str]:
    blob = " ".join(t or "" for t in texts).lower()
    return [k for k in _PROLIF_KW if k in blob]


@lru_cache(maxsize=1024)
def _ot_description(disease: str) -> str:
    try:
        from services.disease_ontology import _ot_graphql, resolve_disease
        d = resolve_disease(disease) or {}
        efo = d.get("efo_id") or d.get("disease_id") or ""
        if not efo:
            return ""
        q = "query($e:String!){ disease(efoId:$e){ description } }"
        r = _ot_graphql(q, {"e": efo}) or {}
        return ((r.get("disease") or {}).get("description") or "")
    except Exception:
        return ""


@lru_cache(maxsize=1024)
def _lit_ratio(disease: str) -> float:
    """PubMed proliferation-literature ratio for a disease — the dynamic catch-all for
    histological hyperproliferation (RA synovial hyperplasia) that text/genetics miss.
    ratio = (disease AND proliferation terms) / (disease). RA≈0.05 vs Alzheimer≈0.00006."""
    try:
        from services import http_client
        def _count(term):
            r = http_client.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                                params={"db": "pubmed", "term": term, "retmode": "json",
                                        "retmax": 0}, timeout=12)
            return int((r.json().get("esearchresult", {}) or {}).get("count", 0)) if r and r.ok else 0
        total = _count(f'"{disease}"[tiab]')
        if total < 20:
            return 0.0
        prolif = _count(f'"{disease}"[tiab] AND (hyperplasia OR "fibroblast proliferation" '
                        'OR hyperproliferation OR pannus OR "synovial proliferation" '
                        'OR "smooth muscle proliferation" OR neointimal)')
        return round(prolif / total, 4)
    except Exception as e:
        logger.debug(f"lit ratio {disease} failed: {e}")
        return 0.0


def hyperproliferative(disease: str, definition: str = "", synonyms: str = "",
                       use_literature: bool = True) -> Dict:
    """{is_prolif, confidence, evidence}. Text signal first (fast); PubMed literature
    ratio as the dynamic catch-all for histological cases text misses."""
    out = {"is_prolif": False, "confidence": 0.0, "evidence": []}
    text_hits = _text_signal(disease, definition, synonyms, _ot_description(disease))
    if text_hits:
        out.update(is_prolif=True, confidence=0.7,
                   evidence=[f"text: {', '.join(text_hits[:3])}"])
    if use_literature:
        ratio = _lit_ratio(disease)
        if ratio >= 0.012:                         # RA≈0.05; benign diseases ≈0.0001
            out["is_prolif"] = True
            out["confidence"] = max(out["confidence"], min(0.9, 0.5 + ratio * 6))
            out["evidence"].append(f"literature ratio {ratio}")
    return out


# ── The combined signal ─────────────────────────────────────────────────────────

def proliferation_match(drug_genes: List[str], drug_actions: Optional[Dict[str, str]],
                        disease: str, definition: str = "", synonyms: str = "") -> Dict:
    """Does the drug INHIBIT a POSITIVE proliferation driver in a HYPERPROLIFERATIVE
    disease? Returns {match, score, driver_gene, direction, disease_confidence, note}.
    Gene gate runs first (cheap, local) so the disease literature query only fires for
    genuine proliferation-inhibitor drugs."""
    out = {"match": False, "score": 0.0, "driver_gene": None, "note": ""}
    actions = {k.upper(): (v or "").lower() for k, v in (drug_actions or {}).items()}
    driver = None
    for g in (drug_genes or []):
        role = proliferation_role(g)
        if not (role["is_driver"] and role["direction"] == "positive"):
            continue
        # Direction of the DRUG on it must be inhibitory (arresting the pathogenic
        # proliferation). If we have no action data, allow it (most cell-cycle drugs
        # ARE inhibitors) but at reduced confidence.
        act = actions.get(g.upper(), "")
        if act and not any(t in act for t in _INHIBITORY):
            continue
        driver = (g.upper(), role, bool(act))
        break
    if not driver:
        return out
    dz = hyperproliferative(disease, definition, synonyms)
    if not dz["is_prolif"]:
        return out
    g, role, has_action = driver
    # score = gene-driver confidence × disease confidence × direction certainty; a
    # small action bonus when the inhibitory action is confirmed. Bounded to a prior.
    score = role["confidence"] * dz["confidence"] * (1.0 if has_action else 0.85)
    out.update(match=True, score=round(min(0.75, score), 4), driver_gene=g,
               direction="positive", disease_confidence=dz["confidence"],
               note=(f"{g} is a positive proliferation driver; the disease is "
                     f"hyperproliferative ({'; '.join(dz['evidence'][:2])}) — a "
                     "cell-cycle inhibitor arrests the pathogenic growth. Credited via "
                     "target-driven proliferation, not genetic-association overlap."))
    return out
