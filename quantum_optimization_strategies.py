"""
Quantum & Geometry Optimization
===============================

Two distinct things, honestly labelled:

1. REAL quantum chemistry — GFN2-xTB (semi-empirical QM) via services.qc_engine,
   when the xTB engine is installed. Returns HOMO/LUMO, gap, dipole,
   polarizability, solvation free energy, QM partial charges and conceptual-DFT
   reactivity indices computed from the electronic structure.

2. Molecular-mechanics 3-D geometry optimization (RDKit UFF/MMFF). This is NOT
   quantum — it is a force field — and is reported as such.

If xTB is unavailable, quantum properties fall back to an RDKit descriptor-based
ESTIMATE that is explicitly flagged `is_estimate: True` so the UI never presents a
heuristic as a real calculation.
"""

import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED, rdPartialCharges, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

try:
    from services import qc_engine
    QC_ENGINE_IMPORTED = True
except Exception:
    QC_ENGINE_IMPORTED = False


class MolecularOptimizer:
    """Runs RDKit geometry optimization and (real or estimated) electronic properties."""

    def optimize(self, smiles: str) -> Dict:
        return optimize_molecular_structure(smiles)

    def get_properties(self, smiles: str) -> Dict:
        return calculate_quantum_properties(smiles)


def optimize_molecular_structure(smiles: str, drug_name: str = "") -> Dict:
    """Molecular-mechanics 3-D conformer optimization (RDKit UFF/MMFF force field).

    Force field, not quantum — labelled accordingly via `ff_type`.
    """
    if not RDKIT_AVAILABLE or not smiles:
        return {"optimized": False, "reason": "RDKit unavailable or no SMILES"}

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"optimized": False, "reason": "Invalid SMILES"}

        mol_h = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3()) != 0:
            if AllChem.EmbedMolecule(mol_h, AllChem.ETDG()) != 0:
                return {"optimized": False, "reason": "Could not embed molecule in 3D"}

        ff = AllChem.UFFGetMoleculeForceField(mol_h)
        if ff is None:
            AllChem.MMFFOptimizeMolecule(mol_h)
            return {"optimized": True, "ff_type": "MMFF",
                    "conformers": mol_h.GetNumConformers(),
                    "num_atoms": mol_h.GetNumAtoms()}

        energy_before = ff.CalcEnergy()
        ff.Minimize(maxIts=500)
        energy_after = ff.CalcEnergy()
        return {
            "optimized": True,
            "ff_type": "UFF",
            "energy_before_kcal": round(energy_before, 3),
            "energy_after_kcal": round(energy_after, 3),
            "energy_delta_kcal": round(energy_before - energy_after, 3),
            "conformers": mol_h.GetNumConformers(),
            "num_atoms": mol_h.GetNumAtoms(),
        }
    except Exception as e:
        logger.warning(f"Geometry optimization failed for {drug_name}: {e}")
        return {"optimized": False, "reason": str(e)}


def _estimate_quantum_properties(smiles: str) -> Dict:
    """RDKit descriptor-based ESTIMATE of electronic properties — used only when the
    real QM engine is unavailable. Every value is flagged as an estimate."""
    if not RDKIT_AVAILABLE or not smiles:
        return {}
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        rdPartialCharges.ComputeGasteigerCharges(mol)
        charges = [a.GetDoubleProp("_GasteigerCharge") for a in mol.GetAtoms()]
        charges = [c for c in charges if c == c]
        max_c, min_c = max(charges, default=0.0), min(charges, default=0.0)
        logp, tpsa = Descriptors.MolLogP(mol), Descriptors.TPSA(mol)
        mw = Descriptors.ExactMolWt(mol)
        gap_proxy = round(max(3.0, min(12.0, 4.5 + 0.01 * tpsa - 0.05 * logp + 0.002 * mw)), 3)
        return {
            "is_estimate": True,
            "gap_ev": gap_proxy,
            "dipole_debye": round((max_c - min_c) * 10, 3),
            "max_partial_charge": round(max_c, 4),
            "min_partial_charge": round(min_c, 4),
            "polarizability_au": round(mol.GetNumHeavyAtoms() * 0.54, 2),
            "qed": round(QED.qed(mol), 4),
            "logp": round(logp, 3),
            "mw": round(mw, 2),
            "tpsa": round(tpsa, 2),
        }
    except Exception as e:
        logger.warning(f"Quantum property estimate failed: {e}")
        return {}


def calculate_quantum_properties(smiles: str, name: str = "") -> Dict:
    """Electronic properties — REAL GFN2-xTB QM when available, else a flagged estimate."""
    if QC_ENGINE_IMPORTED and qc_engine.available():
        res = qc_engine.compute_properties(smiles, name=name)
        if res.get("success"):
            res["is_estimate"] = False
            return res
    return _estimate_quantum_properties(smiles)


def run_quantum_optimization(drug_name: str = "", smiles: str = "") -> Dict:
    """Full pipeline: MM geometry optimization + electronic properties (QM or estimate)."""
    opt = optimize_molecular_structure(smiles, drug_name)
    props = calculate_quantum_properties(smiles, name=drug_name)
    is_real = bool(props) and props.get("is_estimate") is False
    engine = props.get("method") if is_real else "RDKit descriptor estimate (xTB not installed)"
    return {
        "status": "complete" if opt.get("optimized") else "partial",
        "drug_name": drug_name,
        "engine": engine,
        "is_real_qc": is_real,
        "quantum_properties": props,
        "optimization": opt,
        "note": ("Electronic properties from GFN2-xTB semi-empirical quantum mechanics."
                 if is_real else
                 "xTB engine not installed — electronic properties are RDKit descriptor "
                 "ESTIMATES, not a real quantum calculation."),
    }


__all__ = [
    "MolecularOptimizer",
    "optimize_molecular_structure",
    "calculate_quantum_properties",
    "run_quantum_optimization",
]
