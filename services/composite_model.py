"""
Learned composite  —  runtime for the supervised evidence-composite model (gap #4).
═══════════════════════════════════════════════════════════════════════════════
Loads data/composite_model.pkl (fitted by services/train_composite.py) and turns the
evidence sub-signals into a single learned score, REPLACING the hand-weighted sum.
Only present when the learned model actually beat the hand weights (the trainer only
saves it in that case), so: model present → use it; absent → keep the hand composite.

The clinical-reality PENALTIES (safety, CCH, coverage, phantom, registry, trial-failure,
harmful-direction) still apply on top — they are gates on the evidence, not evidence
signals, so they sit outside the learned model. DWPC plausibility likewise stays its
own axis (excluded from the features to avoid circularity).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)

_MODEL_FILE = Path(__file__).parent.parent / "data" / "composite_model.pkl"

# Same order as services/train_composite.FEATURES.
FEATURES = ["target", "pathway", "ppi", "clinical", "indication", "regulatory",
            "network", "directional", "proliferation", "signature", "direction"]

_bundle = None
_loaded = False


def _load():
    global _bundle, _loaded
    if _loaded:
        return
    _loaded = True
    if not _MODEL_FILE.exists():
        logger.info("composite model: not fitted — using hand-weighted composite")
        return
    try:
        import joblib
        _bundle = joblib.load(_MODEL_FILE)
        logger.info("composite model: loaded (cv AUC %s)",
                    (_bundle.get("metrics", {}) or {}).get("auc_learned_logistic_compound_disjoint"))
    except Exception as e:
        logger.warning("composite model load failed: %s", e)
        _bundle = None


def available() -> bool:
    _load()
    return _bundle is not None


def learned_composite(features: List[float]) -> Optional[float]:
    """Learned evidence-composite probability (0..1) from the feature vector (in FEATURES
    order), or None if the model isn't fitted. Fail-soft."""
    _load()
    if _bundle is None:
        return None
    try:
        import numpy as np
        x = _bundle["scaler"].transform(np.asarray([features], dtype=float))
        return float(_bundle["model"].predict_proba(x)[0, 1])
    except Exception as e:
        logger.debug("learned_composite failed: %s", e)
        return None
