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


def is_label_approved_for(drug: str, disease_name: str) -> bool:
    """True when the drug's FDA label indicates it for this disease — i.e. every
    significant disease token appears in the label's indications text. Conservative
    (needs ALL tokens) to avoid false exclusions; fail-soft to False."""
    stems = _disease_stems(disease_name)
    if not stems:
        return False
    text = fda_label_indication_text(drug)
    if not text:
        return False
    # A plural-stemmed token like "arthriti" must still hit the label's "arthritis";
    # substring containment handles the stem→surface-form direction.
    return all(s in text for s in stems)
