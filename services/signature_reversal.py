"""
Signature reversal (CMap / LINCS)  —  the missing orthogonal paradigm (gap #3).
═══════════════════════════════════════════════════════════════════════════════
The most-validated non-target repurposing paradigm (Sirota 2011; CMap/LINCS L1000):
a drug is a candidate for a disease when its transcriptional perturbation REVERSES the
disease's expression signature — it down-regulates what the disease up-regulates and
vice-versa — regardless of whether the drug's target is a disease-association gene.

Our previous signature_engine.py was a PROXY ("assume disease genes are elevated").
This replaces it with REAL differential-expression signatures from CREEDS
(single_drug_perturbations + disease_signatures — GEO-derived up/down gene sets):

  reversal = |drug_down ∩ disease_up| + |drug_up ∩ disease_down|   (drug undoes disease)
  concord  = |drug_up   ∩ disease_up| + |drug_down ∩ disease_down| (drug MIMICS disease)
  connectivity = (reversal − concord) / normaliser   ∈ [−1, 1]

Positive → reversing (a repurposing signal); strongly negative → mimicking (the drug
would push the disease the wrong way — a contraindication-like signal). Fail-soft: no
signature for the drug or disease → covered=False, no score invented.

Coverage: 875 drugs / 828 diseases (crowd-curated, real). Orthogonal to target/KG.
"""
from __future__ import annotations

import json
import logging
import pickle
import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)

_DIR = Path(__file__).parent.parent / "data" / "signatures"
_DRUGS_RAW = _DIR / "creeds_drugs.json"
_DIS_RAW = _DIR / "creeds_diseases.json"
_STORE = _DIR / "signature_store.pkl"          # compact {name -> (up,down)} built once


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (s or "").lower()).strip()


def _genes(lst) -> set:
    """CREEDS gene lists are either ['SYMBOL', ...] or [['SYMBOL', weight], ...]."""
    out = set()
    for g in (lst or []):
        sym = g[0] if isinstance(g, (list, tuple)) and g else g
        if isinstance(sym, str) and sym:
            out.add(sym.upper())                # upper = human/mouse ortholog by symbol
    return out


def _consensus(sigs) -> tuple:
    """Aggregate a name's multiple signatures: a gene is UP if it's up in more of them
    than down (and vice-versa); drop ambiguous genes."""
    up_votes, dn_votes = {}, {}
    for up, dn in sigs:
        for g in up:
            up_votes[g] = up_votes.get(g, 0) + 1
        for g in dn:
            dn_votes[g] = dn_votes.get(g, 0) + 1
    up = {g for g in up_votes if up_votes[g] > dn_votes.get(g, 0)}
    dn = {g for g in dn_votes if dn_votes[g] > up_votes.get(g, 0)}
    return frozenset(up), frozenset(dn)


def _build_store():
    """Parse the two raw CREEDS files into a compact {drug/disease name -> (up,down)}."""
    drug_sigs, dis_sigs = {}, {}
    if _DRUGS_RAW.exists():
        for s in json.loads(_DRUGS_RAW.read_text(encoding="utf-8")):
            for key in filter(None, [s.get("drug_name"), s.get("drugbank_id")]):
                drug_sigs.setdefault(_norm(key), []).append(
                    (_genes(s.get("up_genes")), _genes(s.get("down_genes"))))
    if _DIS_RAW.exists():
        for s in json.loads(_DIS_RAW.read_text(encoding="utf-8")):
            for key in filter(None, [s.get("disease_name"), s.get("do_id"), s.get("umls_cui")]):
                dis_sigs.setdefault(_norm(key), []).append(
                    (_genes(s.get("up_genes")), _genes(s.get("down_genes"))))
    drugs = {k: _consensus(v) for k, v in drug_sigs.items()}
    diseases = {k: _consensus(v) for k, v in dis_sigs.items()}
    _STORE.parent.mkdir(parents=True, exist_ok=True)
    with open(_STORE, "wb") as f:
        pickle.dump({"drugs": drugs, "diseases": diseases}, f)
    logger.info("signature store: %d drugs, %d diseases", len(drugs), len(diseases))
    return {"drugs": drugs, "diseases": diseases}


@lru_cache(maxsize=1)
def _store() -> dict:
    if _STORE.exists():
        try:
            with open(_STORE, "rb") as f:
                return pickle.load(f)
        except Exception as e:
            logger.debug(f"signature store load failed: {e}")
    try:
        return _build_store()
    except Exception as e:
        logger.warning(f"signature store build failed: {e}")
        return {"drugs": {}, "diseases": {}}


def _lookup(name: str, table: dict) -> Optional[tuple]:
    n = _norm(name)
    if n in table:
        return table[n]
    # loose containment for drug/disease name variants
    toks = set(n.split())
    for k, v in table.items():
        kt = set(k.split())
        if toks and (toks <= kt or kt <= toks) and abs(len(n) - len(k)) < 10:
            return v
    return None


def reversal_score(drug: str, disease: str) -> Dict:
    """CMap connectivity for a (drug, disease) pair from real expression signatures.
    {covered, score 0..1 (reversal strength), connectivity −1..1, direction, ...}."""
    out = {"covered": False, "score": 0.0, "connectivity": 0.0, "direction": "none"}
    st = _store()
    dsig = _lookup(drug, st.get("drugs", {}))
    zsig = _lookup(disease, st.get("diseases", {}))
    if not dsig or not zsig:
        return out
    d_up, d_dn = dsig
    z_up, z_dn = zsig
    if not (d_up or d_dn) or not (z_up or z_dn):
        return out
    reversal = len(d_dn & z_up) + len(d_up & z_dn)
    concord = len(d_up & z_up) + len(d_dn & z_dn)
    denom = max(1, reversal + concord)
    conn = (reversal - concord) / denom            # −1..1
    out["covered"] = True
    out["connectivity"] = round(conn, 4)
    out["reversal_hits"] = reversal
    out["concordant_hits"] = concord
    if conn > 0.1:
        out["direction"] = "reversing"
        out["score"] = round(min(1.0, conn), 4)     # repurposing signal
    elif conn < -0.1:
        out["direction"] = "mimicking"              # drug amplifies the disease signature
        out["flag"] = ("Transcriptional signature: the drug MIMICS (does not reverse) "
                       "the disease's expression signature.")
    else:
        out["direction"] = "neutral"
    out["note"] = (f"CMap connectivity {conn:+.2f} from real expression signatures "
                   f"({reversal} reversing vs {concord} concordant gene hits, CREEDS).")
    return out
