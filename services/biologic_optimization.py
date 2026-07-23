"""
Biologics optimization — the antibody/protein counterpart to the SMILES workshop.
═══════════════════════════════════════════════════════════════════════════════
Antibodies (mepolizumab…) are optimized in SEQUENCE + DEVELOPABILITY space, not by
editing atoms. This module returns the modality-appropriate optimization framework:
the developability axes a biologics team actually engineers, each with the design
levers and how it is assessed. When a variable-region SEQUENCE is provided we also
scan the CDRs for the standard chemical-liability motifs (deamidation, isomerization,
oxidation, N-glycosylation, free Cys) — the concrete sequence-level 'edits'.

Honest scope: full affinity-maturation / structure prediction needs specialized
tools + the antibody sequence; this delivers the framework + sequence-liability scan,
clearly labelled, instead of pretending an antibody is a small molecule.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional

# Developability axes an antibody engineering team optimizes (Jain et al. 2017 panel).
_AXES = [
    {"axis": "Affinity / potency",
     "what": "Binding strength to the target epitope (here: IL-5 for mepolizumab).",
     "levers": "Affinity maturation of the CDRs (especially CDR-H3); avoid trading specificity for affinity.",
     "assessed": "SPR/BLI KD, cell-based potency."},
    {"axis": "Specificity",
     "what": "Freedom from off-target / polyspecific binding.",
     "levers": "Reduce positively-charged and hydrophobic CDR patches that drive polyreactivity.",
     "assessed": "PSR / baculovirus-particle (BVP) polyspecificity assays."},
    {"axis": "Thermostability",
     "what": "Resistance to unfolding (shelf-life, manufacturability).",
     "levers": "Framework/CDR stabilizing mutations; strengthen the VH/VL interface; germline the framework.",
     "assessed": "DSF/DSC melting temperature (Tm), Tonset."},
    {"axis": "Solubility",
     "what": "High-concentration solubility for subcutaneous dosing.",
     "levers": "Lower surface hydrophobicity; move the pI away from formulation pH ~6; reduce charge patches.",
     "assessed": "PEG precipitation, cross-interaction chromatography (CIC)."},
    {"axis": "Aggregation",
     "what": "Propensity to form aggregates (immunogenicity + loss of potency).",
     "levers": "Remove aggregation-prone regions (APRs) and exposed hydrophobic patches; IgG1 vs IgG4 format choice.",
     "assessed": "SEC monomer %, accelerated-stability, AC-SINS."},
    {"axis": "Viscosity",
     "what": "Low viscosity at high concentration (needed for auto-injector SC delivery).",
     "levers": "Reduce charge asymmetry (dipole) and hydrophobic patches on the Fv surface.",
     "assessed": "Viscosity at 150 mg/mL; charge-patch / hydrophobic-patch modeling."},
    {"axis": "Immunogenicity",
     "what": "Anti-drug-antibody (ADA) risk.",
     "levers": "Humanize the framework; remove/mask MHC-II T-cell epitopes; germline where possible.",
     "assessed": "In-silico MHC-II epitope prediction, MAPPs / T-cell assays."},
]

# CDR chemical-liability motifs (scanned when a sequence is supplied).
_LIABILITIES = [
    ("N-glycosylation site", r"N[^P][ST]", "Sequon N-x-S/T — non-CDR-desired glycosylation; mutate to remove."),
    ("Deamidation (Asn)", r"N[GSTN]", "NG/NS/NT/NN deamidates → charge variant; substitute the Asn or the +1 residue."),
    ("Isomerization (Asp)", r"D[GSTD]", "DG/DS isomerizes → potency loss; substitute the Asp."),
    ("Oxidation (Met/Trp)", r"[MW]", "Surface Met/Trp oxidizes; consider Leu/Phe where in a CDR contact."),
    ("Free cysteine", r"C", "Unpaired Cys → covalent aggregation/scrambling; check pairing."),
    ("Fragmentation (Asp-Pro)", r"DP", "Acid-labile Asp-Pro bond; fragmentation risk."),
]


# ── Real sequence-derived developability metrics (computed, not hardcoded) ────
# Kyte-Doolittle hydropathy, side-chain pKa, monoisotopic-ish residue MW.
_KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
       "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
       "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
_PKA = {"D": 3.65, "E": 4.25, "C": 8.3, "Y": 10.07, "H": 6.0, "K": 10.53, "R": 12.48}
_POS, _NEG = ("K", "R", "H"), ("D", "E", "C", "Y")
_MW = {"A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14, "Q": 128.13,
       "E": 129.12, "G": 57.05, "H": 137.14, "I": 113.16, "L": 113.16, "K": 128.17,
       "M": 131.19, "F": 147.18, "P": 97.12, "S": 87.08, "T": 101.10, "W": 186.21,
       "Y": 163.18, "V": 99.13}


def _charge_at_pH(seq: str, pH: float) -> float:
    pos = 1.0 / (1.0 + 10 ** (pH - 8.0))          # N-terminus
    for aa in ("K", "R", "H"):
        pos += seq.count(aa) / (1.0 + 10 ** (pH - _PKA[aa]))
    neg = 1.0 / (1.0 + 10 ** (3.1 - pH))          # C-terminus
    for aa in ("D", "E", "C", "Y"):
        neg += seq.count(aa) / (1.0 + 10 ** (_PKA[aa] - pH))
    return pos - neg


def _pI(seq: str) -> float:
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if _charge_at_pH(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def _sequence_metrics(seq: str) -> Dict:
    seq = re.sub(r"[^A-Za-z]", "", seq or "").upper()
    if len(seq) < 8:
        return {}
    n = len(seq)
    gravy = round(sum(_KD.get(a, 0) for a in seq) / n, 3)
    aromatic = round(sum(seq.count(a) for a in "FWY") / n * 100, 1)
    pos = sum(seq.count(a) for a in _POS)
    neg = sum(seq.count(a) for a in _NEG if a != "C")
    mw = round(sum(_MW.get(a, 110) for a in seq) + 18.0, 1)   # +water
    return {
        "length": n,
        "mw_da": mw,
        "pI": _pI(seq),
        "net_charge_pH7": round(_charge_at_pH(seq, 7.4), 1),
        "gravy_hydropathy": gravy,           # >0 = hydrophobic (aggregation/viscosity risk)
        "aromatic_pct": aromatic,
        "pos_residues": pos,
        "neg_residues": neg,
    }


def _metric_flags(m: Dict) -> List[Dict]:
    """Turn computed metrics into developability read-outs (green/amber/red)."""
    if not m:
        return []
    out = []
    pI = m.get("pI")
    if pI is not None:
        # a pI near a common SC formulation pH (~6) risks low solubility there
        lvl = "amber" if 6.5 <= pI <= 8.0 else "green"
        out.append({"metric": "Theoretical pI", "value": pI, "level": lvl,
                    "note": "pI 6.5–8 can hurt solubility at physiological pH; move it away from formulation pH." if lvl == "amber" else "pI is clear of the neutral solubility trough."})
    g = m.get("gravy_hydropathy")
    if g is not None:
        lvl = "red" if g > 0.2 else "amber" if g > -0.2 else "green"
        out.append({"metric": "Hydrophobicity (GRAVY)", "value": g, "level": lvl,
                    "note": "High surface hydrophobicity → aggregation + viscosity risk." if lvl != "green" else "Low hydrophobicity — favourable for solubility."})
    ch = m.get("net_charge_pH7")
    if ch is not None:
        lvl = "amber" if abs(ch) < 1 else "green"
        out.append({"metric": "Net charge @ pH 7.4", "value": ch, "level": lvl,
                    "note": "Near-neutral net charge lowers colloidal stability." if lvl == "amber" else "Net charge supports colloidal stability."})
    ar = m.get("aromatic_pct")
    if ar is not None:
        lvl = "amber" if ar > 12 else "green"
        out.append({"metric": "Aromatic content", "value": str(ar) + "%", "level": lvl,
                    "note": "High Phe/Trp/Tyr content raises hydrophobic-patch/aggregation risk." if lvl == "amber" else "Aromatic content within a typical range."})
    return out


def _scan_sequence(seq: str) -> List[Dict]:
    seq = re.sub(r"[^A-Za-z]", "", seq or "").upper()
    if not seq or len(seq) < 8:
        return []
    hits = []
    for name, pat, note in _LIABILITIES:
        for m in re.finditer(pat, seq):
            hits.append({"liability": name, "motif": m.group(0),
                         "position": m.start() + 1, "note": note})
    # de-dupe common single-residue liabilities to avoid noise (keep first few)
    seen, out = set(), []
    for h in hits:
        k = (h["liability"], h["motif"], h["position"])
        if k in seen:
            continue
        seen.add(k)
        out.append(h)
    return out[:20]


def optimize(name: str, modality: str = "antibody", sequence: str = "") -> Dict:
    """Return the biologics optimization framework for an asset.

    {modality, headline, axes:[...], liabilities:[...], sequence_provided, note}
    """
    is_ab = modality in ("antibody", "fusion_protein")
    liabilities = _scan_sequence(sequence) if sequence else []
    metrics = _sequence_metrics(sequence) if sequence else {}
    metric_flags = _metric_flags(metrics)
    if is_ab:
        headline = (f"{name} is an antibody — optimized in sequence + developability space, "
                    "not by atom-level SMILES edits. Generate variants around the CDR/framework "
                    "regions, then score the developability axes below.")
    else:
        headline = (f"{name} is a {modality.replace('_', ' ')} biologic — optimized by "
                    "sequence/formulation design, not SMILES chemistry.")
    return {
        "modality": modality,
        "headline": headline,
        "axes": _AXES if is_ab else _AXES[2:],   # non-antibody proteins: skip CDR-specific axes
        "liabilities": liabilities,
        "metrics": metrics,
        "metric_flags": metric_flags,
        "sequence_provided": bool(sequence),
        "note": ("Provide the variable-region (VH/VL) sequence to scan the CDRs for chemical "
                 "liabilities (deamidation, isomerization, oxidation, glycosylation, free Cys) — "
                 "the concrete sequence-level edits. Full affinity maturation / structure "
                 "prediction requires dedicated antibody tooling."),
    }


if __name__ == "__main__":
    import json
    r = optimize("mepolizumab", "antibody")
    print(r["headline"])
    for a in r["axes"]:
        print(" -", a["axis"])
