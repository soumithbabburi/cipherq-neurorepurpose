"""
Disk-backed cache (shared)
==========================

A tiny, dependency-free key→value cache that persists under data/cache/<namespace>/.
Replaces the scattered in-memory dicts (protein PDBs, docking results, gene info,
QC properties) that were lost on every restart and duplicated across modules.

Usage:
    cache = DiskCache("protein_pdb", ttl=30 * 86400)   # 30-day TTL
    pdb = cache.get(target)
    if pdb is None:
        pdb = fetch(...)
        cache.set(target, pdb)

Values must be JSON-serialisable (strings, numbers, lists, dicts). Big text blobs
such as PDB files are fine — they are stored as JSON string values.
"""
import hashlib
import json
import logging
import time
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

_CACHE_ROOT = Path(__file__).parent.parent / "data" / "cache"


class DiskCache:
    def __init__(self, namespace: str, ttl: Optional[float] = None):
        """ttl in seconds; None means entries never expire."""
        self.namespace = namespace
        self.ttl = ttl
        self.dir = _CACHE_ROOT / namespace
        try:
            self.dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logger.debug(f"DiskCache mkdir failed for {namespace}: {e}")

    def _path(self, key: str) -> Path:
        h = hashlib.sha1(key.encode("utf-8")).hexdigest()
        return self.dir / f"{h}.json"

    def get(self, key: str) -> Optional[Any]:
        p = self._path(key)
        if not p.exists():
            return None
        try:
            obj = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            return None
        if self.ttl is not None and (time.time() - obj.get("_ts", 0)) > self.ttl:
            return None
        return obj.get("v")

    def set(self, key: str, value: Any) -> None:
        p = self._path(key)
        try:
            tmp = p.with_suffix(".tmp")
            tmp.write_text(json.dumps({"_ts": time.time(), "k": key, "v": value}),
                           encoding="utf-8")
            tmp.replace(p)  # atomic on the same filesystem
        except Exception as e:
            logger.debug(f"DiskCache write failed for {self.namespace}/{key}: {e}")

    def get_or_set(self, key: str, producer) -> Any:
        """Return cached value, or call producer() to compute+store+return it."""
        v = self.get(key)
        if v is not None:
            return v
        v = producer()
        if v is not None:
            self.set(key, v)
        return v
