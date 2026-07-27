# P1 Design: Audit Trail and Authentication

**Roadmap reference:** validation/REGULATORY_POSITIONING.md, P1 items 4 (audit
trail) and 5 (authentication and access control). These are the software
prerequisites for any multi-user or pharma-partner deployment and the core of a
21 CFR Part 11 posture.

**Principle:** additive and reversible. Everything here is gated by an
`AUTH_ENABLED` flag that defaults OFF, so the existing single-user local
workflow is unchanged until the operator opts in. Nothing is deleted; the audit
log is append-only.

Last updated: 2026-07-27.

## 1. Audit trail (Part 11 subpart B — audit trails)

**Requirement.** A secure, computer-generated, time-stamped record of operator
actions that create, modify, or act on records, that cannot be altered without
detection.

**Design.**
* Append-only JSONL at `data/audit/audit_log.jsonl` (a runtime record, gitignored).
* Each entry is HASH-CHAINED: it stores `prev_hash` (the previous entry's hash)
  and its own `hash` = SHA-256 over the canonicalized entry. Editing or deleting
  any past entry breaks the chain from that point on, so tampering is detectable
  by `verify()`. This gives tamper-EVIDENCE without a database.
* Every entry carries: `seq`, UTC `ts`, `actor` (session user, or `anonymous`
  when auth is off), `action`, `details`, and the running `code_commit` (from
  provenance.build_stamp), so an action is traceable to the exact code version.
* Fail-soft: a logging failure is itself logged but never breaks the request.

**Explicitly not yet done:** monthly rotation with chain hand-off, and shipping
the log to append-only external storage (WORM). Noted as future hardening.

## 2. Authentication and access control (Part 11 subpart C — access)

**Requirement.** Limit system access to authorized individuals; use unique
credentials; support role-based authority checks.

**Design.**
* Local user store `data/auth/users.json` (gitignored). Passwords stored only as
  Werkzeug PBKDF2/scrypt hashes, never plaintext.
* Three roles: `admin` (manage users), `analyst` (run screens), `viewer`
  (read-only). Enforced by `login_required` and `role_required` decorators.
* Server-side Flask session. `AUTH_ENABLED` (env, default `false`). When off,
  every route behaves exactly as today and the audit actor is `anonymous`.
* Bootstrap: when auth is enabled and no users exist, a single admin is created
  from `AUTH_ADMIN_USER` + `AUTH_ADMIN_PASSWORD`. There is NO default password;
  if those env vars are absent the app logs a clear error and leaves auth
  unconfigured rather than inventing a weak credential.

**Secret key hardening (prerequisite).** Session security depends on
`app.secret_key`. Today it falls back to a hardcoded string, which is fine for a
local single user but unsafe once sessions carry identity. When auth is enabled
the app requires a real `SECRET_KEY` (env) or a persisted random key at
`data/auth/secret.key`; it refuses to run auth on the hardcoded default.

## 3. Threat model (in scope for this lane)

* Session forgery via a known secret key -> fixed by the secret-key hardening.
* Unauthorized access to screens/records -> fixed by login + role gates.
* Silent tampering with the audit record -> detectable via the hash chain.
* Credential theft at rest -> passwords are hashed (never reversible); the user
  store and secret key are gitignored and must sit outside the repo in a real
  deployment.

Out of scope here (organizational / deployment): TLS termination, SSO/IdP
federation, secrets vaulting, and WORM storage — these belong to the hosting
environment and a partner's QMS, not the application code.

## 4. Verification

* `audit_trail.verify()` walks the chain and returns the first broken sequence
  (or None). Unit-tested: append -> verify OK; a mutated entry -> verify fails.
* `auth` password hash round-trip and role checks are unit-tested.
* Both test sets are pure-logic and run in the CI gate.
