"""
CNS-MPO — Central Nervous System Multiparameter Optimization score
==================================================================
Wager et al., ACS Chem. Neurosci. 2010 — the pharma-standard desirability
metric for CNS (brain-penetrant) drug-likeness. Six physicochemical properties,
each mapped to a 0–1 desirability, summed to 0–6. CNS-MPO >= 4 is the accepted
threshold for a molecule with a good chance of crossing the blood-brain barrier.

For a NEURO repurposing platform this is essential: a non-CNS-penetrant molecule
cannot treat a brain disease regardless of how well it matches the target.

Inputs prefer ChEMBL/ChemAxon computed properties (cx_logp, cx_logd, cx_most_bpka,
full_mwt, psa, hbd); falls back to RDKit from SMILES for molecules without them.
"""
from typing import Dict, Optional

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
    _RDKIT = True
except Exception:
    _RDKIT = False


def _mono_dec(x: float, hi: float, lo: float) -> float:
    """Monotonically decreasing desirability: 1 at x<=hi, 0 at x>=lo, linear between."""
    if x is None:
        return 0.0
    if x <= hi:
        return 1.0
    if x >= lo:
        return 0.0
    return (lo - x) / (lo - hi)


def _hump(x: float, a: float, b: float, c: float, d: float) -> float:
    """Hump desirability: 0 below a, ramps a→b to 1, plateau b→c, ramps c→d to 0."""
    if x is None or x <= a or x >= d:
        return 0.0
    if x < b:
        return (x - a) / (b - a)
    if x <= c:
        return 1.0
    return (d - x) / (d - c)


def _rdkit_props(smiles: str) -> Dict:
    if not (_RDKIT and smiles):
        return {}
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {}
    return {
        "mw":   Descriptors.MolWt(m),
        "logp": Crippen.MolLogP(m),
        "tpsa": rdMolDescriptors.CalcTPSA(m),
        "hbd":  rdMolDescriptors.CalcNumHBD(m),
    }


def cns_mpo(props: Optional[Dict] = None, smiles: str = "") -> Dict:
    """Compute CNS-MPO. `props` may carry ChEMBL keys (full_mwt/mw_freebase,
    cx_logp/alogp, cx_logd, psa, hbd, cx_most_bpka). Missing values fall back to
    RDKit (from `smiles`); a missing basic pKa is treated as non-basic (favourable)."""
    props = props or {}
    rd = _rdkit_props(smiles) if smiles else {}

    def pick(*keys, default=None):
        for k in keys:
            v = props.get(k)
            if v not in (None, "", "None"):
                try:
                    return float(v)
                except (TypeError, ValueError):
                    pass
        return default

    mw   = pick("full_mwt", "mw_freebase", "mw", default=rd.get("mw"))
    clogp = pick("cx_logp", "alogp", "logp", default=rd.get("logp"))
    clogd = pick("cx_logd", default=clogp)          # logD7.4; fall back to logP
    tpsa = pick("psa", "tpsa", default=rd.get("tpsa"))
    hbd  = pick("hbd", default=rd.get("hbd"))
    # Most basic pKa: if absent the molecule is non-basic at physiological pH —
    # favourable for CNS penetration, so score the basicity term as ideal.
    bpka = pick("cx_most_bpka", "most_basic_pka", default=0.0)

    comps = {
        "MW":     {"value": mw,    "desirability": _mono_dec(mw,   360, 500) if mw   is not None else None},
        "cLogP":  {"value": clogp, "desirability": _mono_dec(clogp, 3.0, 5.0) if clogp is not None else None},
        "cLogD":  {"value": clogd, "desirability": _mono_dec(clogd, 2.0, 4.0) if clogd is not None else None},
        "TPSA":   {"value": tpsa,  "desirability": _hump(tpsa, 20, 40, 90, 120) if tpsa is not None else None},
        "HBD":    {"value": hbd,   "desirability": _mono_dec(hbd,  0.5, 3.5) if hbd  is not None else None},
        "pKa(b)": {"value": bpka,  "desirability": _mono_dec(bpka, 8.0, 10.0)},
    }
    # Require the core structural properties (a bare pKa default must not yield a score)
    core = ("MW", "cLogP", "TPSA", "HBD")
    if any(comps[k]["desirability"] is None for k in core):
        return {"score": None, "n_props": 0, "components": comps,
                "cns_druggable": None, "verdict": "insufficient data"}
    avail = [c["desirability"] for c in comps.values() if c["desirability"] is not None]

    # Sum of available desirabilities, scaled to the full 6-property scale so the
    # 4/6 threshold stays meaningful when 1–2 inputs are missing.
    score = round(sum(avail) * (6.0 / len(avail)), 2)
    if score >= 4.5:
        verdict = "High — likely brain-penetrant"
    elif score >= 4.0:
        verdict = "Good — probable CNS penetration"
    elif score >= 3.0:
        verdict = "Borderline — CNS penetration uncertain"
    else:
        verdict = "Low — unlikely to cross the BBB"
    return {
        "score": score, "max": 6, "n_props": len(avail),
        "components": comps, "cns_druggable": score >= 4.0, "verdict": verdict,
    }


if __name__ == "__main__":
    import json
    # Donepezil (CNS drug, should score high) vs erlotinib (oncology, lower CNS)
    print("donepezil:", json.dumps(cns_mpo(smiles="O=C1CC2(CCN(Cc3ccccc3)CC2)Cc2cc(OC)c(OC)cc21"), default=str)[:300])
    print("erlotinib:", json.dumps(cns_mpo(smiles="C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"), default=str)[:300])
