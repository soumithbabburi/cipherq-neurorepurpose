"""
Audit trail (21 CFR Part 11 subpart B) — hash-chained, append-only, rotating.
═══════════════════════════════════════════════════════════════════════════════
A secure, time-stamped, computer-generated record of operator actions that
cannot be altered without detection. Each entry stores the previous entry's hash
and its own SHA-256 over the canonicalized entry, so editing or deleting any
past entry breaks the chain from that point on (tamper-EVIDENCE). Every entry is
stamped with the running code commit so an action traces to the exact version.

The chain is CONTINUOUS across rotated segments: when the active log grows past
_MAX_BYTES it is rotated to a zero-padded, read-only segment file, and the next
entry keeps chaining from the last hash (tracked in a small state pointer), so
verify() validates the entire history across every segment. Read-only segments
are a local tamper-resistance gesture; true WORM/append-only external storage is
a deployment concern (see validation/P1_DESIGN.md).

Fail-soft: a logging failure is recorded to the app log but NEVER breaks the
caller. Storage: data/audit/ (runtime record; gitignored).
"""
from __future__ import annotations

import os
import json
import hashlib
import logging
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DIR = _ROOT / "data" / "audit"
_LOG_FILE = _AUDIT_DIR / "audit_log.jsonl"        # active segment
_GENESIS = "0" * 64                                # prev_hash of the very first entry
_MAX_BYTES = 5_000_000                             # rotate the active segment past this size
_lock = threading.Lock()


# ── paths derived from _AUDIT_DIR at call time (so tests can patch _AUDIT_DIR) ──

def _state_file() -> Path:
    return _AUDIT_DIR / "audit_state.json"


def _segments() -> List[Path]:
    """Rotated segments (in order) followed by the active file, if present."""
    rotated = sorted(_AUDIT_DIR.glob("audit_log.*.jsonl"))
    return rotated + ([_LOG_FILE] if _LOG_FILE.exists() else [])


def _canonical(entry: Dict) -> str:
    payload = {k: entry[k] for k in entry if k != "hash"}
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)


def _hash(entry: Dict) -> str:
    return hashlib.sha256(_canonical(entry).encode("utf-8")).hexdigest()


def _code_commit() -> Optional[str]:
    try:
        from services.provenance import _code_commit as cc
        return cc()
    except Exception:
        return None


def _last_active_entry() -> Optional[Dict]:
    """Last entry in the active file (used only to migrate a missing state pointer)."""
    if not _LOG_FILE.exists():
        return None
    last = None
    with _LOG_FILE.open("r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                last = line.strip()
    try:
        return json.loads(last) if last else None
    except Exception:
        return None


def _read_state() -> Dict:
    """{'last_seq', 'last_hash'} — survives rotation. Falls back to the active file's
    last entry (migration), then to genesis."""
    sf = _state_file()
    if sf.exists():
        try:
            st = json.loads(sf.read_text(encoding="utf-8"))
            if "last_seq" in st and "last_hash" in st:
                return st
        except Exception:
            pass
    prev = _last_active_entry()
    if prev:
        return {"last_seq": prev.get("seq", 0), "last_hash": prev.get("hash", _GENESIS)}
    return {"last_seq": 0, "last_hash": _GENESIS}


def _write_state(seq: int, h: str) -> None:
    _state_file().write_text(json.dumps({"last_seq": seq, "last_hash": h}), encoding="utf-8")


def _rotate_if_needed() -> None:
    """Roll the active segment to a read-only, zero-padded file when it grows past
    _MAX_BYTES. The chain continues because prev_hash comes from the state pointer,
    not from the (now fresh) active file."""
    try:
        if _LOG_FILE.exists() and _LOG_FILE.stat().st_size >= _MAX_BYTES:
            idx = len(sorted(_AUDIT_DIR.glob("audit_log.*.jsonl"))) + 1
            dest = _AUDIT_DIR / f"audit_log.{idx:06d}.jsonl"
            _LOG_FILE.rename(dest)
            try:
                os.chmod(dest, 0o444)      # read-only: local tamper-resistance
            except Exception:
                pass
    except Exception as e:
        logger.error("audit_trail rotation failed: %s", e)


def record(action: str, actor: Optional[str] = None,
           details: Optional[Dict] = None) -> Optional[Dict]:
    """Append one tamper-evident audit entry. Returns the entry, or None if the
    write failed (the failure is logged; the caller is never interrupted)."""
    try:
        with _lock:
            _AUDIT_DIR.mkdir(parents=True, exist_ok=True)
            _rotate_if_needed()
            st = _read_state()
            entry = {
                "seq": st["last_seq"] + 1,
                "ts": datetime.now(timezone.utc).isoformat(),
                "actor": actor or "anonymous",
                "action": str(action),
                "details": details or {},
                "code_commit": _code_commit(),
                "prev_hash": st["last_hash"],
            }
            entry["hash"] = _hash(entry)
            with _LOG_FILE.open("a", encoding="utf-8") as f:
                f.write(json.dumps(entry, default=str) + "\n")
            _write_state(entry["seq"], entry["hash"])
            return entry
    except Exception as e:
        logger.error("audit_trail.record failed for action=%s: %s", action, e)
        return None


def verify() -> Tuple[bool, Optional[int]]:
    """Walk the full chain across all segments in order. Returns (ok,
    first_broken_seq). ok=True on an empty/absent log."""
    prev_hash = _GENESIS
    expected_seq = 1
    try:
        for seg in _segments():
            with seg.open("r", encoding="utf-8") as f:
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
