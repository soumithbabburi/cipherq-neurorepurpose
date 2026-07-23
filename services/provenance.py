"""
Provenance + Live-Data Freshness
════════════════════════════════════════════════════════════════════════════════
Every intelligence platform is only as trustworthy as its sources — and pharma
buyers will not act on a number they cannot trace or date. This layer gives every
data point a LINEAGE (which authoritative source) and a DATA-AGE & INTEGRITY tag
(how fresh, how authoritative), so a live ClinicalTrials.gov pull, a 2023 ChEMBL
snapshot, and a static 2016-lineage knowledge graph are NEVER shown as equally
current.

Design principles (same honesty ethos as the scoring audit):
  • Snapshot dates are GROUNDED in reality — ChEMBL's release date, and the on-disk
    mtime of each local snapshot — never fabricated.
  • LIVE API sources are dated at fetch time (real-time).
  • authority: primary (regulator/registry) > secondary (curated aggregator) >
    derived (model / inferred), so a genetic-association model is never mistaken for
    an FDA approval.
"""
from __future__ import annotations

import os
import logging
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)
_ROOT = Path(__file__).resolve().parent.parent


def _mtime(rel: str) -> Optional[str]:
    """ISO date of a local snapshot file (its real as-of), or None if absent."""
    p = _ROOT / rel
    try:
        if p.exists():
            return date.fromtimestamp(p.stat().st_mtime).isoformat()
    except Exception:
        pass
    return None


def _chembl_release() -> str:
    try:
        from services.neuro_db_service import get_data_footprint
        return ((get_data_footprint() or {}).get("chembl") or {}).get("release") or "2023-05-31"
    except Exception:
        return "2023-05-31"


# ── Source registry ──────────────────────────────────────────────────────────
# live=True  → real-time API (dated at fetch); snapshot_of → local file whose mtime
# is the honest as-of. cadence = how often the UPSTREAM refreshes.
def _registry() -> Dict[str, Dict]:
    return {
        "clinicaltrials": {
            "name": "ClinicalTrials.gov", "kind": "Clinical trials", "authority": "primary",
            "live": True, "cadence": "daily", "url": "https://clinicaltrials.gov/api/v2/studies",
            "integrity": "Registry-verified (trial sponsors, NIH)",
            "feeds": "Trials of a drug in an indication, phase, status, outcomes"},
        "drugs_fda_label": {
            "name": "FDA Drugs@FDA / Labels (openFDA)", "kind": "Regulatory", "authority": "primary",
            "live": True, "cadence": "rolling", "url": "https://api.fda.gov/drug/label.json",
            "integrity": "FDA-approved labelling",
            "feeds": "Approved indications (corrects stale trial-phase → true approval status)"},
        "openfda_faers": {
            "name": "FDA FAERS (openFDA)", "kind": "Pharmacovigilance", "authority": "primary",
            "live": True, "cadence": "quarterly", "url": "https://api.fda.gov/drug/event.json",
            "integrity": "Spontaneous reports — signal, not causation",
            "feeds": "Serious adverse-event signals for the safety / appropriateness gate"},
        "open_targets": {
            "name": "Open Targets Platform", "kind": "Target–disease evidence", "authority": "secondary",
            "live": True, "cadence": "quarterly", "url": "https://api.platform.opentargets.org",
            "integrity": "Curated aggregation (genetics, expression, literature)",
            "feeds": "Target→disease association scores, known-drug indications"},
        "chembl33": {
            "name": "ChEMBL 33", "kind": "Bioactivity", "authority": "secondary",
            "live": False, "cadence": "~yearly", "snapshot": _chembl_release(),
            "url": "https://www.ebi.ac.uk/chembl/", "integrity": "Curated bioactivity DB (local snapshot)",
            "feeds": "Drug identities, targets, measured bioactivities, drug_indication phases"},
        "mondo": {
            "name": "MONDO Disease Ontology", "kind": "Ontology", "authority": "secondary",
            "live": False, "cadence": "monthly", "snapshot": _mtime("data/disease_reference/mondo-simple.json"),
            "url": "https://mondo.monarchinitiative.org/", "integrity": "Community-curated ontology (local snapshot)",
            "feeds": "Disease identities, subtype hierarchy, repurposing-value universe"},
        "drkg": {
            # dated by DATA vintage (DRKG ~2020, Hetionet 2016 lineage), NOT the local
            # file's rebuild time — a recompiled embedding of old data is still old data.
            "name": "DRKG knowledge graph (neurology subset)", "kind": "Knowledge graph", "authority": "derived",
            "live": False, "cadence": "static", "snapshot": "2020-11-01",
            "url": "https://github.com/gnn4dr/DRKG", "integrity": "Model-derived paths — plausibility, not proof; neuro-only coverage; 2016–2020 data vintage",
            "feeds": "Directional mechanism paths + KG plausibility (fail-soft outside coverage)"},
        "orange_book": {
            "name": "FDA Orange Book", "kind": "Regulatory", "authority": "primary",
            "live": False, "cadence": "monthly", "snapshot": _mtime("data/orange_book"),
            "url": "https://www.accessdata.fda.gov/scripts/cder/ob/", "integrity": "FDA exclusivity/patent listings (local snapshot)",
            "feeds": "Small-molecule exclusivity + acquisition signal"},
        "purple_book": {
            "name": "FDA Purple Book", "kind": "Regulatory", "authority": "primary",
            "live": False, "cadence": "monthly", "snapshot": _mtime("data/purple_book"),
            "url": "https://purplebooksearch.fda.gov/", "integrity": "FDA biologic licensure (local snapshot)",
            "feeds": "Biologic reference products for the Biosimilar Radar"},
    }


_FRESH_BANDS = [(0, "fresh", "Live"), (45, "current", "Current"),
                (365, "aging", "Aging"), (10**6, "stale", "Stale")]


def freshness(source_key: str, fetched_at: Optional[str] = None) -> Dict:
    """Data-age classification for a source. Live sources fetched now are 'fresh';
    a snapshot is aged from its real as-of date."""
    reg = _registry().get(source_key, {})
    today = date.today()
    if reg.get("live"):
        as_of = (fetched_at or datetime.now(timezone.utc).isoformat())[:10]
        age = 0
    else:
        as_of = reg.get("snapshot") or None
        if not as_of:
            return {"class": "unknown", "label": "Unknown", "age_days": None, "as_of": None}
        try:
            age = (today - date.fromisoformat(as_of)).days
        except Exception:
            age = None
    cls, label = "unknown", "Unknown"
    if age is not None:
        for thresh, c, lab in _FRESH_BANDS:
            if age <= thresh:
                cls, label = c, lab
                break
    return {"class": cls, "label": label, "age_days": age, "as_of": as_of,
            "live": bool(reg.get("live")), "cadence": reg.get("cadence")}


_AUTH_WEIGHT = {"primary": 1.0, "secondary": 0.75, "derived": 0.5}
_FRESH_WEIGHT = {"fresh": 1.0, "current": 0.9, "aging": 0.65, "stale": 0.4, "unknown": 0.5}


def integrity_tag(source_key: str, fetched_at: Optional[str] = None) -> Dict:
    """Combined DATA-AGE & INTEGRITY indicator for one source: a 0–100 score plus a
    short human tag ('Live · registry-verified', 'Snapshot 2023 · aging', …)."""
    reg = _registry().get(source_key, {})
    fr = freshness(source_key, fetched_at)
    auth = reg.get("authority", "secondary")
    score = round(100 * _AUTH_WEIGHT.get(auth, 0.7) * _FRESH_WEIGHT.get(fr["class"], 0.6))
    if fr.get("live"):
        age_txt = "Live"
    elif fr.get("as_of"):
        age_txt = "Snapshot " + fr["as_of"]
    else:
        age_txt = "Unknown age"
    return {"score": score, "authority": auth, "freshness": fr,
            "tag": f"{age_txt} · {fr['label'].lower()} · {auth}",
            "flag": (fr["class"] in ("stale",)) or (auth == "derived")}


def stamp(source_key: str, fetched_at: Optional[str] = None, **extra) -> Dict:
    """Full provenance record to attach to a datapoint / output field."""
    reg = _registry().get(source_key, {})
    rec = {"source_key": source_key, "source": reg.get("name", source_key),
           "kind": reg.get("kind"), "authority": reg.get("authority"),
           "url": reg.get("url"), "integrity_note": reg.get("integrity"),
           "provenance": integrity_tag(source_key, fetched_at)}
    rec.update(extra)
    return rec


def lineage(source_keys: List[str], fetched_at: Optional[str] = None) -> Dict:
    """Aggregate lineage for an OUTPUT built from several sources: the ordered source
    records + an overall Data-Age & Integrity score (weakest-link aware)."""
    keys = [k for k in source_keys if k in _registry()]
    recs = [stamp(k, fetched_at) for k in keys]
    if not recs:
        return {"sources": [], "overall": None}
    scores = [r["provenance"]["score"] for r in recs]
    overall = round(0.5 * (sum(scores) / len(scores)) + 0.5 * min(scores))  # weakest-link aware
    flags = [r["source"] for r in recs if r["provenance"]["flag"]]
    return {"sources": recs, "overall_score": overall,
            "overall_label": _score_label(overall), "flagged": flags,
            "as_of": (fetched_at or datetime.now(timezone.utc).isoformat())[:10]}


def _score_label(s: int) -> str:
    return "High" if s >= 80 else "Good" if s >= 65 else "Moderate" if s >= 45 else "Low"


def registry_snapshot() -> Dict:
    """All sources with their current freshness + integrity — for the /provenance view."""
    reg = _registry()
    out = []
    for k, v in reg.items():
        it = integrity_tag(k)
        out.append({"key": k, "name": v["name"], "kind": v["kind"], "authority": v["authority"],
                    "url": v["url"], "cadence": v["cadence"], "integrity_note": v["integrity"],
                    "feeds": v.get("feeds", ""), "freshness": it["freshness"],
                    "score": it["score"], "tag": it["tag"], "flag": it["flag"]})
    out.sort(key=lambda r: (-{"primary": 2, "secondary": 1, "derived": 0}.get(r["authority"], 0), -r["score"]))
    live = sum(1 for r in out if r["freshness"].get("live"))
    return {"sources": out, "n": len(out), "live": live,
            "generated_at": datetime.now(timezone.utc).isoformat()}
