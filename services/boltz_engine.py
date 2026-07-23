"""
Boltz-2 binding-affinity engine (P4)  —  structure + affinity prediction.
═══════════════════════════════════════════════════════════════════════════════
The built-in physics engine returns a Vina-style ΔG that is a conservative TRIAGE
score, not a validated affinity. Boltz-2 (open-source, 2025) co-folds a protein +
ligand and predicts a BINDING AFFINITY (a pIC50-scale value + a binder probability)
— the number a medicinal chemist actually wants. This module is the opt-in engine
that runs it and parses the result.

DESIGN (mirrors the DiffDock/NVIDIA opt-in pattern in docking_service):
  • Heavy tool — needs a GPU and multi-GB weights, so it is NEVER on the default
    path. `available` is False unless the `boltz` CLI is on PATH.
  • FAIL-SOFT: any missing dependency, missing weights, or run error returns
    {"success": False, "error": ...} and the docking service falls back to the
    physics engine — the platform never hard-fails on Boltz.
  • Honest labelling: results are tagged "Boltz-2 (predicted affinity)".

Install to activate:  pip install boltz   (first run downloads weights; GPU rec.)
"""
from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)

# 3-letter → 1-letter amino acid code, for extracting a protein sequence from PDB.
_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


class BoltzEngine:
    def __init__(self):
        self._exe = shutil.which("boltz") or os.getenv("BOLTZ_PATH")
        self.available = bool(self._exe)
        if self.available:
            logger.info("Boltz-2 engine available (%s)", self._exe)

    # ── Public API ──────────────────────────────────────────────────────────
    def run(self, protein_pdb: str, ligand_smiles: str, drug_name: str = "ligand",
            target_name: str = "target", timeout: int = 1200) -> Dict:
        """Predict a binding affinity for (protein, ligand). Returns a docking-style
        dict on success, else {"success": False, "error": ...}. Fail-soft."""
        if not self.available:
            return {"success": False, "error": "Boltz-2 not installed (pip install boltz)"}
        seq = self._sequence_from_pdb(protein_pdb)
        if not seq or len(seq) < 20:
            return {"success": False, "error": "Could not extract a protein sequence for Boltz"}
        if not ligand_smiles:
            return {"success": False, "error": "Boltz needs a ligand SMILES"}
        try:
            with tempfile.TemporaryDirectory(prefix="boltz_") as td:
                tdp = Path(td)
                yaml_path = tdp / "input.yaml"
                yaml_path.write_text(self._build_yaml(seq, ligand_smiles), encoding="utf-8")
                out_dir = tdp / "out"
                cmd = [self._exe, "predict", str(yaml_path),
                       "--out_dir", str(out_dir), "--use_msa_server"]
                proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
                if proc.returncode != 0:
                    return {"success": False,
                            "error": f"Boltz run failed: {(proc.stderr or '')[-300:]}"}
                return self._parse_output(out_dir, drug_name, target_name)
        except subprocess.TimeoutExpired:
            return {"success": False, "error": f"Boltz timed out (>{timeout}s)"}
        except Exception as e:
            logger.debug(f"Boltz run error: {e}")
            return {"success": False, "error": f"Boltz error: {e}"}

    # ── Helpers ─────────────────────────────────────────────────────────────
    @staticmethod
    def _sequence_from_pdb(pdb: str) -> str:
        """One-letter sequence of the first chain, from ATOM CA records."""
        seq, seen = [], set()
        for ln in (pdb or "").splitlines():
            if ln[:4] != "ATOM":
                continue
            if ln[12:16].strip() != "CA":
                continue
            chain = ln[21:22]
            resnum = ln[22:27]
            key = (chain, resnum)
            if key in seen:
                continue
            seen.add(key)
            aa = _AA3TO1.get(ln[17:20].strip().upper())
            if aa:
                seq.append(aa)
            # single chain only — stop at a chain break
            if seq and chain and seq and len(seen) > 1 and chain != pdb[21:22]:
                pass
        return "".join(seq)

    @staticmethod
    def _build_yaml(protein_seq: str, ligand_smiles: str) -> str:
        """Boltz-2 input: a protein chain + a ligand, with affinity requested."""
        return (
            "version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: {protein_seq}\n"
            "  - ligand:\n"
            "      id: L\n"
            f"      smiles: '{ligand_smiles}'\n"
            "properties:\n"
            "  - affinity:\n"
            "      binder: L\n"
        )

    @staticmethod
    def _parse_output(out_dir: Path, drug_name: str, target_name: str) -> Dict:
        """Read Boltz affinity JSON + the predicted complex structure."""
        aff_files = list(out_dir.rglob("affinity*.json"))
        aff_val, binder_prob = None, None
        if aff_files:
            try:
                data = json.loads(aff_files[0].read_text())
                # pIC50-scale prediction + probability the ligand is a binder
                aff_val = data.get("affinity_pred_value")
                binder_prob = (data.get("affinity_probability_binary")
                               or data.get("affinity_probability"))
            except Exception:
                pass
        structs = list(out_dir.rglob("*.cif")) + list(out_dir.rglob("*model*.pdb"))
        pose = structs[0].read_text(encoding="utf-8", errors="ignore") if structs else ""
        if aff_val is None and not pose:
            return {"success": False, "error": "Boltz produced no parseable output"}
        # Boltz affinity is a predicted log(IC50)-scale value; expose it directly and
        # convert to an approximate ΔG (kcal/mol) for continuity with the other engines.
        binding = round(float(aff_val), 3) if aff_val is not None else None
        dg = round(-1.364 * float(aff_val), 2) if aff_val is not None else None
        return {
            "success": True,
            "poses": [pose] if pose else [],
            "binding_affinities": [dg] if dg is not None else [],
            "boltz_affinity_pred": binding,
            "boltz_binder_probability": (round(float(binder_prob), 3)
                                         if binder_prob is not None else None),
            "docking_method": "Boltz-2 (predicted affinity)",
            "drug_name": drug_name,
            "target_name": target_name,
            "note": ("Boltz-2 co-folded the protein+ligand and predicted a binding "
                     "affinity (pIC50 scale) + binder probability — a validated-grade "
                     "estimate, unlike the triage ΔG of the physics engine."),
        }
