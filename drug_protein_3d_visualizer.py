"""
Drug-Protein 3D Visualizer using py3Dmol + stmol.
"""

import logging
from typing import Optional

logger = logging.getLogger(__name__)

try:
    import py3Dmol
    PY3DMOL_AVAILABLE = True
except ImportError:
    PY3DMOL_AVAILABLE = False

try:
    from stmol import showmol
    STMOL_AVAILABLE = True
except ImportError:
    STMOL_AVAILABLE = False

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


def _smiles_to_molblock(smiles: str) -> Optional[str]:
    """Convert SMILES to MOL block using RDKit."""
    if not RDKIT_AVAILABLE or not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
        return Chem.MolToMolBlock(mol)
    except Exception as e:
        logger.debug(f"MOL block generation failed: {e}")
        return None


class DrugProtein3DVisualizer:
    """Renders drug molecule and optional protein using py3Dmol + stmol."""

    def __init__(self):
        self.available = PY3DMOL_AVAILABLE and STMOL_AVAILABLE

    def visualize(
        self,
        drug_name: str,
        smiles: str,
        protein_pdb: Optional[str] = None,
        width: int = 800,
        height: int = 450,
    ) -> Optional[object]:
        """
        Build a py3Dmol view of the drug (and optionally the protein).
        Returns the viewer object or None.
        """
        if not PY3DMOL_AVAILABLE:
            return None

        view = py3Dmol.view(width=width, height=height)

        molblock = _smiles_to_molblock(smiles)
        if molblock:
            view.addModel(molblock, "sdf")
            view.setStyle({"model": 0}, {"stick": {"colorscheme": "cyanCarbon", "radius": 0.15}})
            view.addSurface(
                py3Dmol.VDW,
                {"opacity": 0.35, "color": "cyan"},
                {"model": 0},
            )

        if protein_pdb:
            view.addModel(protein_pdb, "pdb")
            view.setStyle({"model": 1}, {"cartoon": {"color": "spectrum", "opacity": 0.8}})

        view.setBackgroundColor("#0a0e1a")
        view.zoomTo()
        return view

    def render_visualization(
        self,
        drug_name: str,
        target_protein: str = "",
        smiles: str = "",
        pdb_content: Optional[str] = None,
        width: int = 800,
        height: int = 450,
    ) -> None:
        """Render the 3D visualization inside a Streamlit context."""
        import streamlit as st

        if not self.available:
            st.info(f"3D visualization requires py3Dmol and stmol. Install with: `pip install py3Dmol stmol`")
            if smiles:
                st.code(f"SMILES: {smiles}")
            return

        view = self.visualize(drug_name, smiles, pdb_content, width=width, height=height)
        if view is None:
            st.warning("Could not generate 3D structure — SMILES may be invalid or RDKit unavailable.")
            return

        showmol(view, height=height, width=width)

    def render(self, drug_data=None, protein_data=None, **kwargs) -> None:
        """Convenience wrapper accepting dicts."""
        import streamlit as st

        drug_name = (drug_data or {}).get("name", "Unknown")
        smiles = (drug_data or {}).get("smiles", "")
        pdb = (protein_data or {}).get("pdb_content")
        target = (protein_data or {}).get("name", "")
        self.render_visualization(drug_name, target, smiles, pdb)

    def render_drug_protein_complex(self, drug_smiles: str, protein_pdb: str, **kwargs) -> None:
        import streamlit as st
        self.render_visualization("Compound", "", drug_smiles, protein_pdb)

    def show_docking_results(self, docking_results, **kwargs) -> None:
        import streamlit as st
        if isinstance(docking_results, dict):
            smiles = docking_results.get("smiles", "")
            name = docking_results.get("drug_name", "")
            if smiles:
                self.render_visualization(name, "", smiles)
            else:
                st.json(docking_results)


visualizer = DrugProtein3DVisualizer()

__all__ = ["DrugProtein3DVisualizer", "visualizer"]
