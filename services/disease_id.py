"""
Disease-ID reconciliation layer.
═══════════════════════════════════════════════════════════════════════════════
The platform pulls disease identities from sources that DO NOT share an id space
or a name convention:
  • Open Targets  → MONDO / EFO ids, canonical lowercase names ("rheumatoid arthritis")
  • ChEMBL        → EFO ids, MeSH inverted headings ("Arthritis, Rheumatoid")
  • openFDA / MeSH→ UMLS CUIs, MeSH headings
  • RepoDB        → UMLS CUIs, free-text indication names

So the SAME disease has different ids and different strings across sources. This
module is the single reconciliation point: it normalises a disease name to a
canonical, order-independent, plural-stemmed token key, and answers whether two
names denote the same disease — independent of ontology or word order. Every
module that matches diseases across sources (approved-indication exclusion, graph
bridging, candidate dedup) should route through here rather than string-compare.

`canonical_key("Arthritis, Rheumatoid") == canonical_key("rheumatoid arthritis")`.
"""
from __future__ import annotations

import re
from typing import FrozenSet

# Non-discriminating disease-name tokens dropped before matching.
_STOP = {
    "disease", "diseases", "disorder", "disorders", "syndrome", "the", "of",
    "and", "with", "chronic", "acute", "type", "primary", "secondary",
    "idiopathic", "nos", "unspecified",
}


def _stem(t: str) -> str:
    """Light plural stemmer (applied to both sides, so consistency is what matters):
    'nephropathies'→'nephropathy', 'neoplasms'→'neoplasm', 'arthritis'→'arthriti'."""
    if len(t) > 4:
        if t.endswith("ies"):
            return t[:-3] + "y"
        if t.endswith("s") and not t.endswith("ss"):
            return t[:-1]
    return t


def canonical_key(name: str) -> FrozenSet:
    """Order-independent, stopword-stripped, plural-stemmed significant-token set
    for a disease name. The platform-wide canonical disease identity for matching."""
    return frozenset(
        _stem(t) for t in re.split(r"[^a-z0-9]+", (name or "").lower())
        if len(t) > 2 and t not in _STOP)


def same_disease(name_a: str, name_b: str) -> bool:
    """True when two disease names denote the same disease: equal canonical token
    sets, or one a subset of the other sharing ≥2 significant tokens (so 'rheumatoid
    arthritis' ≡ 'Arthritis, Rheumatoid' but ≠ 'psoriatic arthritis', and a bare
    'arthritis' does not swallow a specific disease)."""
    a, b = canonical_key(name_a), canonical_key(name_b)
    if not a or not b:
        return False
    if a == b:
        return True
    return (a <= b or b <= a) and len(a & b) >= 2
