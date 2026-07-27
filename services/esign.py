"""
Electronic signatures (21 CFR Part 11 subpart C).
═══════════════════════════════════════════════════════════════════════════════
A Part 11 electronic signature on a record carries three things: the signer's
printed name, the date/time, and the MEANING of the signature (e.g. "reviewed",
"approved"). Signatures here are:

  * executed with two components — the session identity plus the account
    password re-entered at signing time (Part 11 11.200(a)(1));
  * linked to the record via a record_ref and recorded into the tamper-evident,
    hash-chained audit trail, so a signature cannot be excised or transferred
    without breaking the chain (11.70);
  * unique to one individual (the authenticated user).

This is the signing MECHANISM. Binding it to a specific record type (a candidate
dossier, a screen result) is done by the caller passing that record's ref.
Gated by AUTH_ENABLED via auth.verify_credentials. See validation/P1_DESIGN.md
and validation/REGULATORY_POSITIONING.md (P2).
"""
from __future__ import annotations

import logging
from datetime import datetime, timezone
from typing import Dict, List, Optional

from services import auth, audit_trail

logger = logging.getLogger(__name__)

_ACTION = "esignature"


def sign(record_ref: str, meaning: str, username: str, password: str) -> Optional[Dict]:
    """Execute an electronic signature on `record_ref`. Re-authenticates the signer
    with their password (second component) and records the signature manifest into
    the tamper-evident audit trail. Returns the manifest, or None if authentication
    fails or inputs are missing (no signature is recorded on failure)."""
    if not record_ref or not meaning:
        return None
    u = auth.verify_credentials(username, password)
    if not u:
        # record the FAILED signing attempt for the audit trail, but do not sign
        try:
            audit_trail.record("esignature_failed", actor=(username or "").strip().lower(),
                               details={"record_ref": str(record_ref)})
        except Exception:
            pass
        return None
    manifest = {
        "record_ref": str(record_ref),
        "signer": u["username"],           # printed name of the signer
        "signer_role": u["role"],
        "meaning": str(meaning),           # meaning of the signature
        "signed_at": datetime.now(timezone.utc).isoformat(),   # date/time
    }
    audit_trail.record(_ACTION, actor=u["username"], details=manifest)
    return manifest


def signatures_for(record_ref: str) -> List[Dict]:
    """All signatures recorded against a record, oldest first (from the audit trail)."""
    ref = str(record_ref)
    out = []
    for e in audit_trail.iter_entries():
        if e.get("action") == _ACTION and (e.get("details") or {}).get("record_ref") == ref:
            d = dict(e["details"])
            d["audit_seq"] = e.get("seq")      # link back to the chained audit entry
            out.append(d)
    return out
