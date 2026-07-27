"""
Audit trail (21 CFR Part 11 subpart B) — hash-chained, append-only.
═══════════════════════════════════════════════════════════════════════════════
A secure, time-stamped, computer-generated record of operator actions that
cannot be altered without detection. Each entry stores the previous entry's hash
and its own SHA-256 over the canonicalized entry, so editing or deleting any
past entry breaks the chain from that point on (tamper-EVIDENCE). Every entry is
stamped with the running code commit so an action traces to the exact version.

Fail-soft: a logging failure is recorded to the app log but NEVER breaks the
caller. See validation/P1_DESIGN.md.

Storage: data/audit/audit_log.jsonl (runtime record; gitignored).
"""
from __future__ import annotations

import json
import hashlib
import logging
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Tuple

logger = logging.getLogger(__name__)

_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DIR = _ROOT / "data" / "audit"
_LOG_FILE = _AUDIT_DIR / "audit_log.jsonl"
_GENESIS = "0" * 64                       # prev_hash of the first entry
_lock = threading.Lock()


def _canonical(entry: Dict) -> str:
    """Deterministic serialization of the entry MINUS its own hash, for hashing."""
    payload = {k: entry[k] for k in entry if k != "hash"}
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)


def _hash(entry: Dict) -> str:
    return hashlib.sha256(_canonical(entry).encode("utf-8")).hexdigest()


def _last_line() -> Optional[Dict]:
    """The most recent entry (for seq + prev_hash), or None on an empty log."""
    if not _LOG_FILE.exists():
        return None
    last = None
    with _LOG_FILE.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                last = line
    if not last:
        return None
    try:
        return json.loads(last)
    except Exception:
        return None


def _code_commit() -> Optional[str]:
    """Running code version, from the provenance stamp. Never fatal."""
    try:
        from services.provenance import _code_commit as cc
        return cc()
    except Exception:
        return None


def record(action: str, actor: Optional[str] = None,
           details: Optional[Dict] = None) -> Optional[Dict]:
    """Append one tamper-evident audit entry. Returns the entry, or None if the
    write failed (the failure is logged; the caller is never interrupted)."""
    try:
        with _lock:
            _AUDIT_DIR.mkdir(parents=True, exist_ok=True)
            prev = _last_line()
            seq = (prev.get("seq", 0) + 1) if prev else 1
            prev_hash = prev.get("hash", _GENESIS) if prev else _GENESIS
            entry = {
                "seq": seq,
                "ts": datetime.now(timezone.utc).isoformat(),
                "actor": actor or "anonymous",
                "action": str(action),
                "details": details or {},
                "code_commit": _code_commit(),
                "prev_hash": prev_hash,
            }
            entry["hash"] = _hash(entry)
            with _LOG_FILE.open("a", encoding="utf-8") as f:
                f.write(json.dumps(entry, default=str) + "\n")
            return entry
    except Exception as e:
        logger.error("audit_trail.record failed for action=%s: %s", action, e)
        return None


def verify() -> Tuple[bool, Optional[int]]:
    """Walk the chain. Returns (ok, first_broken_seq). ok=True on an empty/absent
    log. first_broken_seq is the seq where the hash or linkage first fails."""
    if not _LOG_FILE.exists():
        return True, None
    prev_hash = _GENESIS
    expected_seq = 1
    try:
        with _LOG_FILE.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                entry = json.loads(line)
                if entry.get("seq") != expected_seq:
                    return False, entry.get("seq")
                if entry.get("prev_hash") != prev_hash:
                    return False, entry.get("seq")
                if _hash(entry) != entry.get("hash"):
                    return False, entry.get("seq")
                prev_hash = entry["hash"]
                expected_seq += 1
    except Exception as e:
        logger.error("audit_trail.verify failed: %s", e)
        return False, expected_seq
    return True, None


def log_file() -> Path:
    return _LOG_FILE
