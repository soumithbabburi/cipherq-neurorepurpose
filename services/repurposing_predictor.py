"""
Repurposing plausibility predictor  —  the validated number.
═══════════════════════════════════════════════════════════════════════════════
Serves P(treats) from the DWPC/metapath logistic model (services/
train_repurposing_predictor.py), the ONE scoring layer with a real, held-out
validation: AUC 0.98 for ranking curated treatments above random pairs, 0.978
compound-disjoint (generalises to unseen drugs).

HONEST SCOPE — this is a MECHANISTIC-PLAUSIBILITY / triage score, NOT a probability
of clinical success. On RepoDB it separates real treatments from random noise
(AUC 0.98) but NOT eventual successes from eventual failures (AUC 0.42): both are
biologically plausible, and efficacy/safety/PK are not in the graph. The
data-derived actionable cutoff (Youden-optimal on the treats-vs-random task) is
P ≥ 0.15 — "plausible enough for a human to review", not "likely to work".

COVERAGE — only the Hetionet subgraph (8,894 DrugBank compounds, 213 DO diseases)
can be scored. Off-coverage pairs return None (fail-soft) and the caller keeps the
evidence score. A precomputed DWPC cache makes lookups O(1) after first load.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, Optional

import numpy as np

logger = logging.getLogger(__name__)

_ROOT = Path(__file__).parent.parent
_MODEL_FILE = _ROOT / "data" / "repurposing_predictor.pkl"
_CACHE_FILE = _ROOT / "data" / "dwpc_cache.npz"

ACTIONABLE_CUTOFF = 0.15          # Youden-optimal P(treats) for treats-vs-random


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (s or "").lower()).strip()


# Tokens that carry no discriminating power for disease/drug name matching.
_STOP = {"disease", "disorder", "syndrome", "the", "of", "and", "s", "type",
         "chronic", "acute", "primary", "secondary", "malignant", "benign"}


def _tokens(s: str) -> frozenset:
    return frozenset(t for t in _norm(s).split() if t and t not in _STOP)


class _Predictor:
    def __init__(self):
        self.ok = False
        self.clf = self.scaler = None
        self.feats = []
        self.stack = None                       # (n_cmp, n_dis, n_feat) float32
        self.cmp_name2i: Dict[str, int] = {}
        self.dis_name2i: Dict[str, int] = {}
        self.dis_names: list = []               # disease index -> Hetionet display name
        self.cmp_names: list = []               # compound index -> Hetionet display name
        self._loaded = False

    def _load(self):
        if self._loaded:
            return
        self._loaded = True
        try:
            import joblib
            b = joblib.load(_MODEL_FILE)
            self.clf, self.scaler, self.feats = b["model"], b["scaler"], b["features"]
        except Exception as e:
            logger.info(f"plausibility predictor: model unavailable ({e})")
            return
        try:
            self._load_or_build_cache()
            self.ok = self.stack is not None
        except Exception as e:
            logger.warning(f"plausibility predictor: DWPC cache failed ({e})")

    def _load_or_build_cache(self):
        if _CACHE_FILE.exists():
            d = np.load(_CACHE_FILE, allow_pickle=True)
            # "dis_names" was added for the KG candidate GENERATOR — rebuild older caches.
            if list(d["features"]) == list(self.feats) and "dis_names" in d.files:
                self.stack = d["stack"]
                self.cmp_name2i = d["cmp_name2i"].item()
                self.dis_name2i = d["dis_name2i"].item()
                self.dis_names = list(d["dis_names"])
                self.cmp_names = list(d["cmp_names"])
                return
        self._build_cache()

    def _build_cache(self):
        from services.metapath_features import get_features_engine
        import psycopg2
        from config import db_params
        eng = get_features_engine(log=lambda *a, **k: None)
        mats = [eng.dwpc[f] for f in self.feats]
        self.stack = np.stack(mats, axis=-1).astype(np.float32)   # (n_cmp,n_dis,n_feat)

        # node id → matrix index, then NAME → index for resolution + index → display name
        cmp_id2i, dis_id2i = eng.cmp_index, eng.dis_index
        n_cmp, n_dis = self.stack.shape[0], self.stack.shape[1]
        self.cmp_names = [""] * n_cmp
        self.dis_names = [""] * n_dis
        conn = psycopg2.connect(**db_params()); cur = conn.cursor()
        cur.execute("SELECT id, name, kind FROM hetionet_nodes "
                    "WHERE kind IN ('Compound','Disease')")
        for nid, nm, kind in cur.fetchall():
            if not nm:
                continue
            key = _norm(nm)
            if kind == "Compound" and nid in cmp_id2i:
                self.cmp_name2i.setdefault(key, cmp_id2i[nid])
                self.cmp_names[cmp_id2i[nid]] = nm
            elif kind == "Disease" and nid in dis_id2i:
                self.dis_name2i.setdefault(key, dis_id2i[nid])
                self.dis_names[dis_id2i[nid]] = nm
        conn.close()
        try:
            np.savez_compressed(_CACHE_FILE, stack=self.stack, features=np.array(self.feats),
                                cmp_name2i=self.cmp_name2i, dis_name2i=self.dis_name2i,
                                dis_names=np.array(self.dis_names, dtype=object),
                                cmp_names=np.array(self.cmp_names, dtype=object))
        except Exception as e:
            logger.debug(f"DWPC cache save skipped: {e}")

    def _resolve(self, name: str, table: Dict[str, int]) -> Optional[int]:
        n = _norm(name)
        if not n:
            return None
        if n in table:                                  # 1. exact normalized match
            return table[n]
        qt = _tokens(name)
        if not qt:
            return None
        # 2. token match: the candidate whose meaningful tokens best overlap the
        #    query — require full containment either way, or high Jaccard, so
        #    "alzheimer disease" ↔ "alzheimer s disease" resolves but unrelated
        #    names do not. Ties broken toward the tightest (fewest extra tokens).
        best = None                                     # (score, extra, idx)
        for k, i in table.items():
            kt = _tokens(k)
            if not kt:
                continue
            inter = len(qt & kt)
            if not inter:
                continue
            contained = qt <= kt or kt <= qt
            jac = inter / len(qt | kt)
            if contained or jac >= 0.6:
                cand = (jac, -(len(qt ^ kt)), i)
                if best is None or cand > best:
                    best = cand
        return best[2] if best else None

    def predict(self, drug_name: str, disease_name: str) -> Optional[Dict]:
        self._load()
        if not self.ok:
            return None
        ci = self._resolve(drug_name, self.cmp_name2i)
        di = self._resolve(disease_name, self.dis_name2i)
        if ci is None or di is None:
            return None
        x = np.log1p(self.stack[ci, di][None, :].astype(float))
        p = float(self.clf.predict_proba(self.scaler.transform(x))[0, 1])
        return {
            "probability": round(p, 4),
            "actionable": p >= ACTIONABLE_CUTOFF,
            "cutoff": ACTIONABLE_CUTOFF,
            "covered": True,
            "basis": "DWPC metapath model (Hetionet); validated AUC 0.98 vs random",
            "scope": ("Mechanistic plausibility / triage - NOT probability of clinical "
                      "success (AUC 0.42 vs real trial failures)."),
        }

    # ── KG-as-candidate-generator (gap #2) ───────────────────────────────────────
    # Instead of only SCORING candidates that OT target-associations already surfaced,
    # the validated DWPC model can GENERATE them: score the drug against EVERY Hetionet
    # disease (one matrix pass) and return the top P(treats). This is the Rephetio/TxGNN
    # inversion — the KG model drives discovery, evidence enriches — and unblocks the
    # narrow OT-bounded universe (e.g. Palbociclib's 25 candidates).
    def top_diseases_for_drug(self, drug_name: str, k: int = 30,
                              min_p: float = ACTIONABLE_CUTOFF) -> List[Dict]:
        self._load()
        if not self.ok:
            return []
        ci = self._resolve(drug_name, self.cmp_name2i)
        if ci is None:
            return []
        X = np.log1p(self.stack[ci].astype(float))          # (n_dis, n_feat)
        p = self.clf.predict_proba(self.scaler.transform(X))[:, 1]
        order = np.argsort(-p)
        out = []
        for i in order[:k * 2]:
            if p[i] < min_p:
                break
            nm = self.dis_names[i] if i < len(self.dis_names) else ""
            if nm:
                out.append({"disease": nm, "probability": round(float(p[i]), 4)})
            if len(out) >= k:
                break
        return out

    def top_drugs_for_disease(self, disease_name: str, k: int = 30,
                              min_p: float = ACTIONABLE_CUTOFF) -> List[Dict]:
        self._load()
        if not self.ok:
            return []
        di = self._resolve(disease_name, self.dis_name2i)
        if di is None:
            return []
        X = np.log1p(self.stack[:, di].astype(float))       # (n_cmp, n_feat)
        p = self.clf.predict_proba(self.scaler.transform(X))[:, 1]
        order = np.argsort(-p)
        out = []
        for i in order[:k * 2]:
            if p[i] < min_p:
                break
            nm = self.cmp_names[i] if i < len(self.cmp_names) else ""
            if nm:
                out.append({"drug": nm, "probability": round(float(p[i]), 4)})
            if len(out) >= k:
                break
        return out


_singleton = _Predictor()


def plausibility(drug_name: str, disease_name: str) -> Optional[Dict]:
    """P(treats) for a (drug, disease) pair, or None off Hetionet coverage. Fail-soft."""
    try:
        return _singleton.predict(drug_name, disease_name)
    except Exception as e:
        logger.debug(f"plausibility predict failed: {e}")
        return None


def top_diseases_for_drug(drug_name: str, k: int = 30) -> List[Dict]:
    """KG-generated candidate INDICATIONS for a drug — the drug's highest-P(treats)
    diseases across the whole Hetionet, ranked. [] off coverage. Fail-soft."""
    try:
        return _singleton.top_diseases_for_drug(drug_name, k)
    except Exception as e:
        logger.debug(f"top_diseases_for_drug failed: {e}")
        return []


def top_drugs_for_disease(disease_name: str, k: int = 30) -> List[Dict]:
    """KG-generated candidate DRUGS for a disease — highest-P(treats) compounds. Fail-soft."""
    try:
        return _singleton.top_drugs_for_disease(disease_name, k)
    except Exception as e:
        logger.debug(f"top_drugs_for_disease failed: {e}")
        return []


if __name__ == "__main__":
    import sys
    d = sys.argv[1] if len(sys.argv) > 1 else "Metformin"
    dis = sys.argv[2] if len(sys.argv) > 2 else "Alzheimer's disease"
    print(f"{d} -> {dis}:", plausibility(d, dis))
