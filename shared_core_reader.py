"""
shared_core_reader.py — read the certified shared-core shelf (CipherQ side).
════════════════════════════════════════════════════════════════════════════
Twin of POZ's reader. The orchestration pipeline (../POZ-main/orchestration)
publishes ChEMBL + OpenFDA + ClinicalTrials + Open Targets to the `shared_core`
schema inside the shared chembl_33 Postgres. This lets CipherQ READ that
certified data instead of re-hitting the same live APIs.

Read-through cache, not a hard dependency: every function is defensive and
returns empty on any error, so the app always falls back to its own live path.

Connection: reuses CipherQ's own config.db_params() but pins the database to
chembl_33 (where shared_core lives), exactly like neuro_db_service does.
"""
_ENGINE = None
DEFAULT_MAX_AGE_DAYS = 30


def _engine():
    global _ENGINE
    if _ENGINE is None:
        from urllib.parse import quote_plus
        from sqlalchemy import create_engine
        import config
        p = config.db_params()
        pw = quote_plus(str(p.get("password", "")))
        auth = f"{p['user']}:{pw}@" if pw else f"{p['user']}@"
        # shared_core lives in chembl_33, not the default neurorepurpose DB.
        _ENGINE = create_engine(
            f"postgresql+psycopg2://{auth}{p['host']}:{p['port']}/chembl_33",
            pool_pre_ping=True,
        )
    return _ENGINE


def get_opentargets_scores(disease_name, max_age_days=DEFAULT_MAX_AGE_DAYS):
    """Certified {drug_name_lower: phase_score} for a disease, or {} if not
    cached / stale / unavailable."""
    if not disease_name or not disease_name.strip():
        return {}
    try:
        from sqlalchemy import text as sa_text
        with _engine().connect() as c:
            rows = c.execute(
                sa_text(
                    "SELECT drug_name, score FROM shared_core.opentargets_drugs "
                    "WHERE lower(disease) = lower(:d) "
                    "AND certified_at > now() - make_interval(days => :age)"
                ),
                {"d": disease_name.strip(), "age": int(max_age_days)},
            ).fetchall()
        return {r[0]: float(r[1]) for r in rows}
    except Exception:
        return {}


def get_trials(drug_name, max_age_days=DEFAULT_MAX_AGE_DAYS):
    """Certified trial rows for a drug (list of dicts), or [] if not cached."""
    if not drug_name or not drug_name.strip():
        return []
    try:
        from sqlalchemy import text as sa_text
        with _engine().connect() as c:
            rows = c.execute(
                sa_text(
                    "SELECT nct_id, title, status, phase, sponsor, conditions "
                    "FROM shared_core.clinical_trials "
                    "WHERE lower(query_drug) = lower(:d) "
                    "AND certified_at > now() - make_interval(days => :age)"
                ),
                {"d": drug_name.strip(), "age": int(max_age_days)},
            ).fetchall()
        return [
            {"nct_id": r[0], "title": r[1], "status": r[2], "phase": r[3],
             "sponsor": r[4], "conditions": r[5]}
            for r in rows
        ]
    except Exception:
        return []


def is_available():
    try:
        from sqlalchemy import text as sa_text
        with _engine().connect() as c:
            c.execute(sa_text("SELECT 1 FROM shared_core.certification LIMIT 1"))
        return True
    except Exception:
        return False
