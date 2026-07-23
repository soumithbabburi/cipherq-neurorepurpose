"""
Train the PoS calibrator  →  data/pos_model.pkl
═══════════════════════════════════════════════════════════════════════════════
Fits a real, supervised model on outcomes already in the platform's database, so
the "AI/ML" in Asset Scoring is an actual fitted classifier with a documented
schema — not the orphaned ml_scoring_model.pkl that nothing loaded.

Label   : did the molecule reach approval (max_phase == 4) for ANY indication?
Features: biological-profile signals the platform computes anyway —
            n_targets        distinct mechanistic targets
            max_pchembl      strongest measured potency (pChEMBL)
            n_indications    breadth of indication evidence
            n_mechanisms     count of annotated mechanisms

NOTE on leakage: the current clinical `phase` is deliberately EXCLUDED from the
features. Because the label is derived from phase (approved = phase >= 4),
including phase would let the model trivially memorise the label (ROC-AUC ≈ 1.0)
— a classic target leak. Phase is handled rigorously by the analytic
phase-transition model in pos_model.py; this calibrator contributes the
ORTHOGONAL biological-approvability signal, which is what makes the 50/50 blend
informative rather than circular.

The output bundle is {"model", "features", "meta"} and is consumed by
services/pos_model.py as an optional calibration blend. If the database is not
reachable, this script explains how to run it and exits cleanly — the analytic
PoS model keeps working without it.

Usage:  python -m services.train_pos_model        (from the repo root)
"""
from __future__ import annotations

import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger("train_pos")

_ROOT = Path(__file__).parent.parent
_OUT = _ROOT / "data" / "pos_model.pkl"

# `phase` is intentionally NOT a feature — it derives the label (target leak).
FEATURES = ["n_targets", "max_pchembl", "n_indications", "n_mechanisms"]

# One query: per-compound features + the approval label, straight from ChEMBL-
# derived tables already loaded for scoring.
_SQL = """
WITH base AS (
    SELECT c.id AS compound_id,
           COALESCE(MAX(i.max_phase), 0)                  AS phase,
           COUNT(DISTINCT m.target_id)                    AS n_targets,
           COALESCE(MAX(a.pchembl_value), 0)              AS max_pchembl,
           COUNT(DISTINCT i.id)                           AS n_indications,
           COUNT(DISTINCT m.id)                           AS n_mechanisms
    FROM compounds c
    LEFT JOIN indications i        ON i.compound_id = c.id
    LEFT JOIN mechanisms  m        ON m.compound_id = c.id
    LEFT JOIN compound_activities a ON a.compound_id = c.id
    GROUP BY c.id
)
SELECT phase, n_targets, max_pchembl, n_indications, n_mechanisms,
       CASE WHEN phase >= 4 THEN 1 ELSE 0 END AS approved
FROM base
WHERE n_targets > 0 OR n_indications > 0;
"""


def _fetch_rows():
    try:
        import psycopg2
        import psycopg2.extras
        from config import db_params
    except Exception as e:
        logger.warning("DB libraries/config unavailable: %s", e)
        return None
    try:
        conn = psycopg2.connect(**db_params())
    except Exception as e:
        logger.warning("Could not connect to the database: %s", e)
        return None
    try:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(_SQL)
            return [dict(r) for r in cur.fetchall()]
    finally:
        conn.close()


def train() -> bool:
    rows = _fetch_rows()
    if not rows:
        logger.info(
            "\nNo training data (database not reachable).\n"
            "The analytic phase-transition PoS model works without this file.\n"
            "To fit the calibrator, start Postgres and run:\n"
            "    python -m services.train_pos_model\n"
        )
        return False

    try:
        import joblib
        import numpy as np
        from sklearn.calibration import CalibratedClassifierCV
        from sklearn.ensemble import GradientBoostingClassifier
        from sklearn.model_selection import cross_val_score
    except Exception as e:
        logger.warning("scikit-learn not installed: %s  (pip install scikit-learn joblib)", e)
        return False

    X = np.array([[float(r[f]) for f in FEATURES] for r in rows], dtype=float)
    y = np.array([int(r["approved"]) for r in rows], dtype=int)
    pos = int(y.sum())
    logger.info("Training rows: %d  (approved=%d, %.1f%%)", len(y), pos, 100 * pos / max(1, len(y)))
    if pos < 20 or (len(y) - pos) < 20:
        logger.warning("Too few examples in one class to train a reliable model — skipping.")
        return False

    base = GradientBoostingClassifier(random_state=0, max_depth=3, n_estimators=150)
    clf = CalibratedClassifierCV(base, method="isotonic", cv=3)
    try:
        auc = float(cross_val_score(base, X, y, cv=5, scoring="roc_auc").mean())
        logger.info("5-fold ROC-AUC: %.3f", auc)
    except Exception:
        auc = None
    clf.fit(X, y)

    # ── Integrity gate ────────────────────────────────────────────────────────
    # An implausibly high AUC on this label means the features are leaking it.
    # Here, annotation richness (n_indications/n_mechanisms/...) is a CONSEQUENCE
    # of a drug being approved and marketed, not information available at the
    # decision point — reverse causation. We refuse to certify such a model as
    # trustworthy, so pos_model.py will ignore it and use the analytic model.
    LEAK, FLOOR = 0.95, 0.55
    trustworthy = auc is not None and FLOOR <= auc <= LEAK
    if not trustworthy:
        logger.warning(
            "\n⚠ Calibrator NOT certified for blending (cv_auc=%s).\n"
            "  AUC>%.2f indicates the label is leaking via reverse causation\n"
            "  (ChEMBL annotation richness reflects approval, it does not predict\n"
            "  it). The analytic phase-transition model is used instead. Saving the\n"
            "  bundle flagged untrustworthy for transparency/audit.",
            None if auc is None else round(auc, 3), LEAK)

    bundle = {
        "model": clf,
        "features": FEATURES,
        "meta": {
            "n_rows": len(y), "n_approved": pos, "cv_auc": auc,
            "trustworthy": trustworthy,
            "label": "max_phase >= 4 (reached approval for any indication)",
            "note": "Blended 50/50 with the analytic phase-transition model in "
                    "services/pos_model.py ONLY when trustworthy is True.",
        },
    }
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, _OUT)
    logger.info("Saved calibrator → %s (trustworthy=%s)", _OUT, trustworthy)
    return True


if __name__ == "__main__":
    train()
