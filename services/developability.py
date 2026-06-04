"""
Area-aware Developability Scorer

Drug-likeness is not one-size-fits-all. A topical dermatology candidate is judged
on stratum-corneum permeation; an eye drop on corneal permeation; a CNS drug on
blood–brain-barrier penetration; an oral systemic drug on Lipinski's Rule-of-Five.

Given a molecule's SMILES and a therapeutic area (or an explicit area filter), this
module picks the right profile from data/developability_profiles.json, computes the
relevant physicochemical descriptors with RDKit, and scores how well the molecule
fits the properties that matter FOR THAT ROUTE.

The profiles/ranges live in the JSON config (editable) — no thresholds are hardcoded
here. RDKit is optional; without it the module degrades to "unavailable" cleanly.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

PROFILES_FILE = Path(__file__).parent.parent / "data" / "developability_profiles.json"

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDKIT_OK = True
except Exception:
    RDKIT_OK = False

_config_cache: Optional[dict] = None


def _config() -> dict:
    global _config_cache
    if _config_cache is None:
        try:
            _config_cache = json.loads(PROFILES_FILE.read_text(encoding="utf-8"))
        except Exception as e:
            logger.warning(f"Could not load developability profiles: {e}")
            _config_cache = {"area_to_profile": {"_default": "oral_systemic"}, "profiles": {}}
    return _config_cache


def is_available() -> bool:
    return RDKIT_OK and bool(_config().get("profiles"))


# ── Profile selection (dynamic from therapeutic area) ───────────────────────

def profile_for_area(area: str = "", therapeutic_areas: Optional[List[str]] = None) -> str:
    """Map an area filter and/or a disease's Open Targets therapeutic areas to a profile key."""
    cfg = _config()
    mapping = cfg.get("area_to_profile", {})
    default = mapping.get("_default", "oral_systemic")

    candidates = []
    if area:
        candidates.append(area.lower())
    for ta in (therapeutic_areas or []):
        candidates.append(ta.lower())

    for cand in candidates:
        for key, profile in mapping.items():
            if key == "_default":
                continue
            if key in cand or cand in key:
                return profile
    return default


# ── Descriptor computation ──────────────────────────────────────────────────

def descriptors(smiles: str) -> Dict[str, float]:
    """Compute the physicochemical descriptors used by the profiles."""
    if not RDKIT_OK or not smiles:
        return {}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    return {
        "mw":   round(Descriptors.MolWt(mol), 1),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "tpsa": round(Descriptors.TPSA(mol), 1),
        "hbd":  rdMolDescriptors.CalcNumHBD(mol),
        "hba":  rdMolDescriptors.CalcNumHBA(mol),
        "rtb":  rdMolDescriptors.CalcNumRotatableBonds(mol),
        "qed":  round(Descriptors.qed(mol), 3),
    }


def _evaluate(value: float, spec: dict) -> tuple:
    """Return (status, status_value) for a property value against its spec."""
    imin = spec.get("ideal_min", float("-inf"))
    imax = spec.get("ideal_max", float("inf"))
    amin = spec.get("acceptable_min", imin if "ideal_min" in spec else float("-inf"))
    amax = spec.get("acceptable_max", imax if "ideal_max" in spec else float("inf"))

    if imin <= value <= imax:
        return ("ideal", 1.0)
    if amin <= value <= amax:
        return ("acceptable", 0.6)
    return ("out", 0.0)


def _target_text(spec: dict) -> str:
    unit = spec.get("unit", "")
    lo = spec.get("ideal_min")
    hi = spec.get("ideal_max")
    if lo is not None and hi is not None:
        return f"{lo}–{hi}{unit}"
    if hi is not None:
        return f"≤ {hi}{unit}"
    if lo is not None:
        return f"≥ {lo}{unit}"
    return "—"


# ── Public scoring API ──────────────────────────────────────────────────────

def score(smiles: str, area: str = "", therapeutic_areas: Optional[List[str]] = None) -> Dict:
    """
    Score a molecule's developability against the area-appropriate profile.

    Returns: {
      available, profile, profile_label, route, rationale,
      score (0-1), level (High/Moderate/Low), properties: [
        {key, label, value, unit, target, status}  # status: ideal|acceptable|out
      ], descriptors: {...}
    }
    """
    cfg = _config()
    profile_key = profile_for_area(area, therapeutic_areas)
    profile = (cfg.get("profiles") or {}).get(profile_key, {})

    base = {
        "available":     False,
        "profile":       profile_key,
        "profile_label": profile.get("label", profile_key),
        "route":         profile.get("route", ""),
        "rationale":     profile.get("rationale", ""),
        "score":         None,
        "level":         None,
        "properties":    [],
        "descriptors":   {},
    }

    if not is_available():
        return base
    desc = descriptors(smiles)
    if not desc:
        return base

    props_spec = profile.get("properties", {})
    rows: List[Dict] = []
    total_w = 0.0
    earned  = 0.0
    for key, spec in props_spec.items():
        if key not in desc:
            continue
        val = desc[key]
        status, sval = _evaluate(val, spec)
        w = spec.get("weight", 1)
        total_w += w
        earned  += w * sval
        rows.append({
            "key":    key,
            "label":  spec.get("label", key),
            "value":  val,
            "unit":   spec.get("unit", ""),
            "target": _target_text(spec),
            "status": status,
        })

    s = round(earned / total_w, 3) if total_w else None
    level = None
    if s is not None:
        level = "High" if s >= 0.75 else "Moderate" if s >= 0.5 else "Low"

    base.update({
        "available":   True,
        "score":       s,
        "level":       level,
        "properties":  rows,
        "descriptors": desc,
    })
    return base
