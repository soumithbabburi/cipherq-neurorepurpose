"""
Quantum Optimization Strategies
Real molecular optimization using RDKit with Streamlit UI rendering.
"""

import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import (
        Descriptors, QED, rdMolDescriptors, AllChem,
        rdPartialCharges, rdForceFieldHelpers,
    )
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


class MolecularOptimizer:
    """Runs RDKit-based molecular optimization and quantum property estimation."""

    def optimize(self, smiles: str) -> Dict:
        return optimize_molecular_structure(smiles)

    def get_properties(self, smiles: str) -> Dict:
        return calculate_quantum_properties(smiles)


def optimize_molecular_structure(smiles: str, drug_name: str = "") -> Dict:
    """
    3D-conformer optimization via UFF force field (RDKit).
    Returns energy delta and conformer count.
    """
    if not RDKIT_AVAILABLE or not smiles:
        return {"optimized": False, "reason": "RDKit unavailable or no SMILES"}

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"optimized": False, "reason": "Invalid SMILES"}

        mol_h = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
        if result != 0:
            # Fallback: random distance geometry
            result = AllChem.EmbedMolecule(mol_h, AllChem.ETDG())

        if result != 0:
            return {"optimized": False, "reason": "Could not embed molecule in 3D"}

        ff = AllChem.UFFGetMoleculeForceField(mol_h)
        if ff is None:
            ff_type = "MMFF"
            AllChem.MMFFOptimizeMolecule(mol_h)
        else:
            ff_type = "UFF"
            energy_before = ff.CalcEnergy()
            ff.Minimize(maxIts=500)
            energy_after = ff.CalcEnergy()
            energy_delta = energy_before - energy_after
            return {
                "optimized": True,
                "ff_type": ff_type,
                "energy_before_kcal": round(energy_before, 3),
                "energy_after_kcal": round(energy_after, 3),
                "energy_delta_kcal": round(energy_delta, 3),
                "conformers": mol_h.GetNumConformers(),
                "num_atoms": mol_h.GetNumAtoms(),
            }

        return {"optimized": True, "ff_type": ff_type, "conformers": mol_h.GetNumConformers()}

    except Exception as e:
        logger.warning(f"Optimization failed for {drug_name}: {e}")
        return {"optimized": False, "reason": str(e)}


def calculate_quantum_properties(smiles: str) -> Dict:
    """
    Estimate quantum-chemical properties using RDKit descriptors.
    Provides HOMO-LUMO gap proxy, partial charge spread, and electronic properties.
    """
    if not RDKIT_AVAILABLE or not smiles:
        return {}

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}

        # Gasteiger partial charges as electron density proxy
        rdPartialCharges.ComputeGasteigerCharges(mol)
        charges = [atom.GetDoubleProp("_GasteigerCharge") for atom in mol.GetAtoms()]
        charges = [c for c in charges if c == c]  # filter NaN

        max_charge = max(charges, default=0.0)
        min_charge = min(charges, default=0.0)
        charge_spread = max_charge - min_charge  # proxy for dipole moment

        # HOMO-LUMO gap proxy: ionisation potential ≈ -E_HOMO estimate from TPSA + LogP
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        mw = Descriptors.ExactMolWt(mol)
        qed = QED.qed(mol)

        # Empirical HOMO-LUMO proxy (eV range typical for drug-like molecules: 4–10 eV)
        homo_lumo_proxy = 4.5 + 0.01 * tpsa - 0.05 * logp + 0.002 * mw
        homo_lumo_proxy = round(max(3.0, min(12.0, homo_lumo_proxy)), 3)

        # Polarisability proxy (Å³): ~0.54 * heavy atoms
        polarisability = round(mol.GetNumHeavyAtoms() * 0.54, 2)

        return {
            "homo_lumo_gap_eV": homo_lumo_proxy,
            "dipole_moment_proxy_D": round(charge_spread * 10, 3),
            "max_partial_charge": round(max_charge, 4),
            "min_partial_charge": round(min_charge, 4),
            "polarisability_A3": polarisability,
            "qed": round(qed, 4),
            "logp": round(logp, 3),
            "mw": round(mw, 2),
            "tpsa": round(tpsa, 2),
        }

    except Exception as e:
        logger.warning(f"Quantum property calculation failed: {e}")
        return {}


def run_quantum_optimization(drug_name: str = "", smiles: str = "") -> Dict:
    """Run full optimization pipeline and return combined results."""
    opt = optimize_molecular_structure(smiles, drug_name)
    props = calculate_quantum_properties(smiles)
    return {
        "status": "complete" if opt.get("optimized") else "partial",
        "drug_name": drug_name,
        "optimization": opt,
        "quantum_properties": props,
    }


def render_quantum_optimization_section(drug_name: str = "", smiles: str = "", drug_data: Dict = None) -> None:
    """Render quantum optimization results in Streamlit."""
    import streamlit as st

    if drug_data and not drug_name:
        drug_name = drug_data.get("name", "")
    if drug_data and not smiles:
        smiles = drug_data.get("smiles", "")

    st.subheader("Quantum & 3D Optimization")

    if not smiles:
        st.info("No SMILES structure available for this compound.")
        return

    with st.spinner("Computing molecular properties…"):
        results = run_quantum_optimization(drug_name, smiles)

    opt = results.get("optimization", {})
    props = results.get("quantum_properties", {})

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**3D Conformer Optimization**")
        if opt.get("optimized"):
            st.success(f"Optimised with {opt.get('ff_type','UFF')} force field")
            if "energy_delta_kcal" in opt:
                st.metric("Energy reduction", f"{opt['energy_delta_kcal']:.2f} kcal/mol")
            st.metric("Conformers generated", opt.get("conformers", 1))
        else:
            st.warning(f"Optimization limited: {opt.get('reason', 'unknown')}")

    with col2:
        st.markdown("**Electronic Properties**")
        if props:
            st.metric("HOMO-LUMO gap (proxy)", f"{props.get('homo_lumo_gap_eV', 0):.2f} eV")
            st.metric("Dipole moment (proxy)", f"{props.get('dipole_moment_proxy_D', 0):.2f} D")
            st.metric("Polarisability", f"{props.get('polarisability_A3', 0):.1f} Å³")
        else:
            st.info("Electronic properties unavailable")

    if props:
        with st.expander("Full quantum property data"):
            st.json(props)


__all__ = [
    "MolecularOptimizer",
    "optimize_molecular_structure",
    "calculate_quantum_properties",
    "run_quantum_optimization",
    "render_quantum_optimization_section",
]
