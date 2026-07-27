"""
Authentication and access control (21 CFR Part 11 subpart C).
═══════════════════════════════════════════════════════════════════════════════
Local user store with Werkzeug-hashed passwords, three roles, and Flask
session-based login. GATED by AUTH_ENABLED (env, default false): when off, every
route behaves exactly as before and the audit actor is 'anonymous'. No default
password is ever invented. See validation/P1_DESIGN.md.

Storage: data/auth/users.json (runtime secret; gitignored).
"""
from __future__ import annotations

import os
import json
import logging
import threading
from datetime import datetime, timezone
from functools import wraps
from pathlib import Path
from typing import Dict, Optional

from werkzeug.security import generate_password_hash, check_password_hash

logger = logging.getLogger(__name__)

_ROOT = Path(__file__).resolve().parent.parent
_AUTH_DIR = _ROOT / "data" / "auth"
_USERS_FILE = _AUTH_DIR / "users.json"
_SECRET_FILE = _AUTH_DIR / "secret.key"
_HARDCODED_DEFAULT_SECRET = "neurorepurpose-cipherq-2026"   # must NOT back real sessions

ROLES = ("admin", "analyst", "viewer")
_ROLE_RANK = {"viewer": 0, "analyst": 1, "admin": 2}
_lock = threading.Lock()


def auth_enabled() -> bool:
    return os.getenv("AUTH_ENABLED", "false").strip().lower() in ("1", "true", "yes", "on")


# ── User store ────────────────────────────────────────────────────────────────

def _load_users() -> Dict[str, Dict]:
    if not _USERS_FILE.exists():
        return {}
    try:
        return json.loads(_USERS_FILE.read_text(encoding="utf-8")) or {}
    except Exception as e:
        logger.error("auth: user store unreadable: %s", e)
        return {}


def _save_users(users: Dict[str, Dict]) -> None:
    _AUTH_DIR.mkdir(parents=True, exist_ok=True)
    _USERS_FILE.write_text(json.dumps(users, indent=2), encoding="utf-8")


def add_user(username: str, password: str, role: str = "viewer") -> bool:
    """Create or replace a user with a hashed password. Returns False on bad input."""
    username = (username or "").strip().lower()
    if not username or not password or role not in ROLES:
        return False
    with _lock:
        users = _load_users()
        users[username] = {
            "pw_hash": generate_password_hash(password),
            "role": role,
            "created": datetime.now(timezone.utc).isoformat(),
        }
        _save_users(users)
    return True


def set_password(username: str, password: str) -> bool:
    username = (username or "").strip().lower()
    with _lock:
        users = _load_users()
        if username not in users or not password:
            return False
        users[username]["pw_hash"] = generate_password_hash(password)
        _save_users(users)
    return True


def verify_credentials(username: str, password: str) -> Optional[Dict]:
    """Return {'username', 'role'} on a correct password, else None."""
    username = (username or "").strip().lower()
    user = _load_users().get(username)
    if not user or not password:
        return None
    if check_password_hash(user.get("pw_hash", ""), password):
        return {"username": username, "role": user.get("role", "viewer")}
    return None


def has_role(user_role: Optional[str], required: str) -> bool:
    """True if user_role is at least as privileged as `required`."""
    if user_role not in _ROLE_RANK or required not in _ROLE_RANK:
        return False
    return _ROLE_RANK[user_role] >= _ROLE_RANK[required]


def list_users() -> list:
    """Users WITHOUT password hashes — for the admin view / audit."""
    return sorted(({"username": u, "role": v.get("role", "viewer"),
                    "created": v.get("created")} for u, v in _load_users().items()),
                  key=lambda r: r["username"])


def required_role_for(method: str, path: str) -> str:
    """Declarative authority policy (Part 11): the minimum role for a request.
    admin for /admin/*, analyst for any state-changing verb, viewer for reads.
    One auditable place instead of decorating dozens of routes."""
    p = path or ""
    if p.startswith("/admin"):
        return "admin"
    if (method or "GET").upper() in ("POST", "PUT", "PATCH", "DELETE"):
        return "analyst"
    return "viewer"


def bootstrap_admin() -> Optional[str]:
    """When auth is enabled and no users exist, create a single admin from
    AUTH_ADMIN_USER + AUTH_ADMIN_PASSWORD. Returns the username, or None if it
    could not (missing env vars) — never invents a default password."""
    if not auth_enabled() or _load_users():
        return None
    u = (os.getenv("AUTH_ADMIN_USER") or "").strip().lower()
    p = os.getenv("AUTH_ADMIN_PASSWORD") or ""
    if not u or not p:
        logger.error("auth: AUTH_ENABLED but no users and AUTH_ADMIN_USER/"
                     "AUTH_ADMIN_PASSWORD not set — auth is unconfigured.")
        return None
    add_user(u, p, "admin")
    logger.info("auth: bootstrapped admin user '%s'", u)
    return u


# ── Secret-key hardening ────────────────────────────────────────────────────────

def resolve_secret_key() -> Optional[str]:
    """A strong session secret: env SECRET_KEY (if not the hardcoded default),
    else a persisted random key at data/auth/secret.key (created once). Returns
    None only if the filesystem is unwritable AND no safe env key is set."""
    env = os.getenv("SECRET_KEY", "")
    if env and env != _HARDCODED_DEFAULT_SECRET:
        return env
    try:
        if _SECRET_FILE.exists():
            k = _SECRET_FILE.read_text(encoding="utf-8").strip()
            if k:
                return k
        _AUTH_DIR.mkdir(parents=True, exist_ok=True)
        k = os.urandom(32).hex()
        _SECRET_FILE.write_text(k, encoding="utf-8")
        logger.info("auth: generated a persistent random SECRET_KEY at %s", _SECRET_FILE)
        return k
    except Exception as e:
        logger.error("auth: could not resolve a strong secret key: %s", e)
        return None


# ── Flask decorators (import flask lazily so this module is testable standalone) ──

def current_user() -> Optional[Dict]:
    try:
        from flask import session
        u = session.get("user")
        return u if isinstance(u, dict) else None
    except Exception:
        return None


def login_required(fn):
    @wraps(fn)
    def wrapper(*a, **k):
        if not auth_enabled():
            return fn(*a, **k)
        from flask import redirect, url_for, request
        if not current_user():
            return redirect(url_for("login", next=request.path))
        return fn(*a, **k)
    return wrapper


def role_required(required: str):
    def deco(fn):
        @wraps(fn)
        def wrapper(*a, **k):
            if not auth_enabled():
                return fn(*a, **k)
            from flask import redirect, url_for, request, abort
            u = current_user()
            if not u:
                return redirect(url_for("login", next=request.path))
            if not has_role(u.get("role"), required):
                abort(403)
            return fn(*a, **k)
        return wrapper
    return deco
