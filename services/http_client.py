"""
Shared HTTP client (retry + backoff + sane timeouts)
====================================================

Every external-API module (Open Targets, ChEMBL, STRING, NCBI, ClinicalTrials,
RCSB, UniProt, AlphaFold) previously rolled its own `requests` call with
inconsistent (or missing) timeouts and no retry. One transient blip then quietly
degraded results. This module centralises that: a pooled session, a bounded
timeout, and exponential-backoff retry on transient failures (timeouts, conn
errors, 429/5xx).

    from services.http_client import get_json, post_json, get_text
    data = get_json(url, params={...}, default={})
"""
import logging
import time
from typing import Any, Optional

import requests

logger = logging.getLogger(__name__)

_RETRY_STATUS = {429, 500, 502, 503, 504}
_session: Optional[requests.Session] = None


def _get_session() -> requests.Session:
    global _session
    if _session is None:
        s = requests.Session()
        s.headers.update({"User-Agent": "RepurposeIQ/1.0"})
        _session = s
    return _session


def request(method: str, url: str, *, timeout: float = 8.0, retries: int = 2,
            backoff: float = 0.6, **kw) -> Optional[requests.Response]:
    """Issue a request, retrying transient failures with exponential backoff.

    Returns the Response (which may carry a non-retryable error status), or None
    if every attempt raised a network exception.
    """
    last_exc: Optional[Exception] = None
    for attempt in range(retries + 1):
        try:
            r = _get_session().request(method, url, timeout=timeout, **kw)
            if r.status_code in _RETRY_STATUS and attempt < retries:
                time.sleep(backoff * (2 ** attempt))
                continue
            return r
        except requests.RequestException as e:
            last_exc = e
            if attempt < retries:
                time.sleep(backoff * (2 ** attempt))
                continue
    if last_exc:
        logger.debug(f"{method} {url} failed after {retries + 1} attempts: {last_exc}")
    return None


def get(url: str, **kw) -> Optional[requests.Response]:
    return request("GET", url, **kw)


def post(url: str, **kw) -> Optional[requests.Response]:
    return request("POST", url, **kw)


def get_json(url: str, *, default: Any = None, **kw) -> Any:
    r = get(url, **kw)
    if r is not None and r.status_code == 200:
        try:
            return r.json()
        except ValueError:
            logger.debug(f"GET {url}: non-JSON response")
    return default


def post_json(url: str, *, default: Any = None, **kw) -> Any:
    r = post(url, **kw)
    if r is not None and r.status_code == 200:
        try:
            return r.json()
        except ValueError:
            logger.debug(f"POST {url}: non-JSON response")
    return default


def get_text(url: str, *, default: str = "", **kw) -> str:
    r = get(url, **kw)
    if r is not None and r.status_code == 200:
        return r.text
    return default
