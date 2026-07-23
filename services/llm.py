"""
Local LLM helper (Ollama / llama3.1:8b).
════════════════════════════════════════════════════════════════════════════════
A thin, fail-soft client for the locally running Ollama server. Used only to phrase
GROUNDED facts into readable prose — never to invent facts. Every caller passes the
real, computed evidence in the prompt and instructs the model to stay within it, so
a model hiccup can only cost us fluency, never correctness (callers keep a
deterministic fallback string).

House style is enforced after generation, not just requested: the platform copy uses
no emojis, no decorative symbols, and no hyphens or dashes, so `clean_style()` strips
them from whatever the model returns.
"""
from __future__ import annotations

import logging
import re
from typing import Optional

logger = logging.getLogger(__name__)

OLLAMA_URL = "http://127.0.0.1:11434/api/generate"
MODEL = "llama3.1:8b"

_available: Optional[bool] = None


def available() -> bool:
    """Is the Ollama server reachable and holding our model? Cached after first check."""
    global _available
    if _available is not None:
        return _available
    try:
        import requests
        r = requests.get("http://127.0.0.1:11434/api/tags", timeout=2)
        _available = bool(r.ok and MODEL.split(":")[0] in r.text)
    except Exception as e:
        logger.debug("Ollama not available: %s", e)
        _available = False
    return _available


# Emoji + pictographic ranges (kept explicit so the intent is auditable).
_EMOJI = re.compile(
    "[\U0001F300-\U0001FAFF\U00002600-\U000027BF\U0001F1E6-\U0001F1FF"
    "\U00002190-\U000021FF\U00002B00-\U00002BFF️‍]"
)
# Dashes and en/em variants → a plain space (house rule: no hyphens or dashes).
_DASH = re.compile(r"\s*[‐-―−\-]+\s*")


def clean_style(text: str) -> str:
    """Enforce house style: no emojis, no markdown/decorative symbols, no dashes."""
    if not text:
        return ""
    t = _EMOJI.sub("", text)
    t = t.replace("**", "").replace("*", "").replace("`", "").replace("#", "")
    t = t.replace("’", "'").replace("‘", "'").replace("“", '"').replace("”", '"')
    t = _DASH.sub(" ", t)                       # hyphen/en/em/minus → space
    t = re.sub(r"[ \t]{2,}", " ", t)
    t = re.sub(r"\n{3,}", "\n\n", t)
    return t.strip()


def generate(prompt: str, system: str = "", *, max_tokens: int = 320,
             temperature: float = 0.2, timeout: int = 45) -> Optional[str]:
    """Generate text from the local model, style-cleaned. Returns None on any failure
    (server down, timeout, bad response) so callers can fall back deterministically."""
    if not available():
        return None
    try:
        import requests
        payload = {
            "model": MODEL,
            "prompt": prompt,
            "stream": False,
            "options": {"temperature": temperature, "num_predict": max_tokens},
        }
        if system:
            payload["system"] = system
        r = requests.post(OLLAMA_URL, json=payload, timeout=timeout)
        if not r.ok:
            return None
        out = (r.json() or {}).get("response", "")
        cleaned = clean_style(out)
        return cleaned or None
    except Exception as e:
        logger.debug("LLM generate failed: %s", e)
        return None
