"""
Approval supplement  —  catch ChEMBL's per-indication phase lag.
═══════════════════════════════════════════════════════════════════════════════
ChEMBL's `max_phase_for_ind` lags real FDA approvals (e.g. Baricitinib was
approved for atopic dermatitis / alopecia areata years before ChEMBL marks those
indications phase 4), so the reverse engine's approved-indication exclusion —
which keys on ChEMBL phase ≥ 4 — lets a drug's OWN newly-approved indication
resurface as a "novel" candidate.

This module supplements that with the authoritative per-indication source: the
current FDA label's INDICATIONS AND USAGE section (openFDA drug/label). NOTE: the
FDA Orange Book (products.txt) lists approved PRODUCTS + dates + patents but NOT
the indication each treats — the indication lives in the label — so the label is
the correct source here.

FAIL-SOFT: no label / network error → empty text → no exclusion added (we never
block a candidate on missing data).
"""
from __future__ import annotations

import logging
import re
from typing import FrozenSet

from services import http_client

logger = logging.getLogger(__name__)

OPENFDA_LABEL = "https://api.fda.gov/drug/label.json"

_text_cache: dict = {}


def fda_label_indication_text(drug: str) -> str:
    """Lowercased INDICATIONS AND USAGE text from the current FDA label for a drug.
    Cached per drug; '' on any failure."""
    d = (drug or "").strip().lower()
    if not d:
        return ""
    if d in _text_cache:
        return _text_cache[d]
    text = ""
    # Try each name field in turn (openFDA treats "+"-joined clauses as AND, so a
    # single combined query requiring brand==generic returns nothing). First field
    # that resolves a label wins.
    for field in ("generic_name", "brand_name", "substance_name"):
        try:
            r = http_client.get(OPENFDA_LABEL, params={
                "search": f'openfda.{field}:"{d}"', "limit": 5}, timeout=12)
            if r is not None and r.ok:
                for rec in r.json().get("results", []):
                    iu = " ".join(rec.get("indications_and_usage", []) or [])
                    if iu.strip():
                        text += " " + iu.lower()
            if text:
                break
        except Exception as e:
            logger.debug(f"FDA label fetch failed for {drug} ({field}): {e}")
    _text_cache[d] = text
    return text


def _disease_stems(disease_name: str) -> FrozenSet:
    """Reuse the reverse engine's normalized disease tokens (order-independent,
    plural-stemmed) so label matching is consistent with the rest of the platform."""
    try:
        from services.reverse_repurposing import _disease_tokens
        return _disease_tokens(disease_name)
    except Exception:
        import re
        return frozenset(t for t in re.split(r"[^a-z0-9]+", (disease_name or "").lower())
                         if len(t) > 2)


def _label_names_disease(text: str, stems: FrozenSet) -> bool:
    """True iff the candidate's significant token set `stems` appears in `text` as a
    SELF-CONTAINED disease phrase — a run of exactly those tokens (interleaved only with
    stopwords) that is NOT directly preceded by another significant token which would make
    it a QUALIFIED variant of a different disease.

    Why (bug 2026-07): the old check `all(tok in text)` let a BROADER disease be marked
    approved by a more-SPECIFIC approved string that merely contains its name — e.g.
    "granulomatosis with polyangiitis" (GPA) matched the label's "EOSINOPHILIC
    granulomatosis with polyangiitis" (EGPA, a different disease). Here the extra qualifier
    "eosinophilic" sits directly before the phrase, so that run is rejected; GPA no longer
    matches. Fully general — no disease names, no qualifier list; the discriminating token
    is whatever significant word the more-specific indication carries and the candidate lacks.
    """
    from services.disease_id import _STOP, _stem
    S = set(stems)
    if not S:
        return False
    raw = [w for w in re.split(r"[^a-z0-9]+", (text or "").lower()) if w]
    toks = [(_stem(w), (len(w) > 2 and w not in _STOP)) for w in raw]   # (stem, is_significant)
    n = len(toks)
    for a in range(n):
        st, is_sig = toks[a]
        if not is_sig or st not in S:
            continue
        # LEFT flank: an adjacent preceding SIGNIFICANT token absent from the candidate means
        # the label's disease extends left into a qualified variant (e.g. eosinophilic|GPA).
        if a > 0 and toks[a - 1][1] and toks[a - 1][0] not in S:
            continue
        seen, ok, b = set(), True, a
        while b < n and len(seen) < len(S):
            bst, bsig = toks[b]
            if bsig:
                if bst in S:
                    seen.add(bst)
                else:                    # a foreign significant token interrupts the phrase
                    ok = False
                    break
            b += 1
        if ok and seen == S:
            return True
    return False


def is_label_approved_for(drug: str, disease_name: str) -> bool:
    """True when the drug's FDA label indicates it for THIS disease. Requires the disease's
    significant tokens to appear as a self-contained phrase in the label's indications text
    (not merely as a substring of a more-specific, qualified indication). Fail-soft to False."""
    stems = _disease_stems(disease_name)
    if not stems:
        return False
    text = fda_label_indication_text(drug)
    if not text:
        return False
    return _label_names_disease(text, stems)
