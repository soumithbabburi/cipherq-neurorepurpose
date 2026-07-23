"""
Vina docking engine (app-side)  —  real AutoDock Vina, native Windows.
═══════════════════════════════════════════════════════════════════════════════
Bridges the Flask app (running in its own venv) to a real AutoDock Vina docking
run. The cheminformatics stack (RDKit/Meeko/pdbfixer/OpenBabel) lives ONLY in the
`drugdisc-agent` conda env, and AutoDock Vina ships as a native Windows CLI binary
(the Python API has no win-64 build), so this engine simply:

  1. writes the receptor PDB + a JSON job to a temp dir,
  2. shells out to vendor/vina_worker.py under the drugdisc-agent python,
  3. reads back the poses (SDF) + real ΔG (kcal/mol).

Unlike the built-in physics engine (a Vina-STYLE empirical approximation), this is
the real thing: proper receptor protonation, Meeko ligand PDBQT, and Vina's own
search + scoring, docked INTO the co-crystal binding site when the structure is
holo. Fail-soft: `available` is False when the env/binary/worker are missing, so
the docking service transparently falls back to the physics engine.
"""
from __future__ import annotations

import json
import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

logger = logging.getLogger(__name__)

_ROOT = Path(__file__).resolve().parent.parent

# Overridable via env; defaults match the repo layout / the drugdisc-agent env
# created alongside the quantum engine's micromamba root at .qc/.
_DRUGDISC_PY = Path(os.getenv(
    "DRUGDISC_PYTHON",
    _ROOT / ".qc" / "envs" / "drugdisc-agent" / "python.exe"))
_VINA_EXE = Path(os.getenv(
    "VINA_EXE", _ROOT / "vendor" / "vina_bin" / "vina.exe"))
_WORKER = _ROOT / "vendor" / "vina_worker.py"


def _clean_receptor(pdb: str) -> str:
    """Strip a receptor PDB down to protein atoms only before docking:
      - drop ALL HETATM (co-crystal ligand, waters, ions) — critically, a holo
        structure's bound ligand SITS IN the pocket we dock into; leaving it there
        blocks the box and can crash the Vina worker (receptor-prep access violation).
        The dock box centre is already fixed from pocket detection, so removing the
        ligand is both the crash fix and the scientifically correct empty-pocket prep.
      - keep only the blank/'A' alternate location so duplicate altloc atoms don't
        confuse protonation.
    """
    out = []
    for ln in pdb.splitlines():
        rec = ln[:6].strip()
        if rec == "ATOM":
            alt = ln[16] if len(ln) > 16 else " "
            if alt not in (" ", "A"):
                continue
            out.append(ln)
        elif rec == "TER":
            out.append(ln)
    out.append("END")
    return "\n".join(out)


class VinaEngine:
    def __init__(self):
        self.available = (_DRUGDISC_PY.exists() and _VINA_EXE.exists()
                          and _WORKER.exists())
        if self.available:
            logger.info("VinaEngine ready — AutoDock Vina 1.2.7 (native)")
        else:
            missing = [str(p) for p in (_DRUGDISC_PY, _VINA_EXE, _WORKER)
                       if not p.exists()]
            logger.info(f"VinaEngine unavailable (missing: {missing}) — "
                        "docking will use the physics engine")

    def run(self, protein_pdb: str, smiles: str, drug_name: str,
            target_name: str, dock_center: Optional[List[float]] = None,
            exhaustiveness: int = 8, n_poses: int = 9,
            timeout: int = 600) -> Dict:
        """Dock `smiles` into `protein_pdb` with AutoDock Vina. Returns the standard
        docking-service dict (poses=SDF strings, binding_affinities=kcal/mol)."""
        if not self.available:
            return {"success": False, "error": "Vina engine unavailable", "poses": []}
        if not protein_pdb or len(protein_pdb) < 200:
            return {"success": False, "error": "No receptor structure for Vina", "poses": []}

        with tempfile.TemporaryDirectory(prefix="vina_job_") as td:
            tdp = Path(td)
            prot_file = tdp / "receptor.pdb"
            prot_file.write_text(_clean_receptor(protein_pdb), encoding="utf-8")
            out_json = tdp / "result.json"
            job = {
                "protein_pdb_file": str(prot_file),
                "smiles": smiles,
                "vina_exe": str(_VINA_EXE),
                "out_json": str(out_json),
                "exhaustiveness": exhaustiveness,
                "n_poses": n_poses,
            }
            if dock_center is not None:
                job["center"] = [float(c) for c in dock_center]
            job_file = tdp / "job.json"
            job_file.write_text(json.dumps(job), encoding="utf-8")

            try:
                proc = subprocess.run(
                    [str(_DRUGDISC_PY), str(_WORKER), str(job_file)],
                    capture_output=True, text=True, timeout=timeout)
            except subprocess.TimeoutExpired:
                return {"success": False, "error": f"Vina timed out after {timeout}s", "poses": []}

            if not out_json.exists():
                return {"success": False,
                        "error": f"Vina worker produced no result "
                                 f"(rc={proc.returncode}): {(proc.stderr or '')[-400:]}",
                        "poses": []}
            res = json.loads(out_json.read_text(encoding="utf-8"))

        if not res.get("success"):
            return {"success": False, "error": res.get("error", "Vina failed"),
                    "poses": [], "trace": res.get("trace", "")}

        affs = res.get("binding_affinities", []) or []
        # Map ΔG → 0–1 confidence, same convention as the physics engine so the UI
        # bars read consistently (≈ −4 kcal/mol → low, ≈ −12 → high).
        confs = [round(float(1.0 / (1.0 + np.exp((a + 7.0) / 1.5))), 3) for a in affs]
        return {
            "success": True,
            "poses": res["poses"],
            "binding_affinities": affs,
            "confidence_scores": confs,
            "docking_method": res.get("docking_method", "AutoDock Vina 1.2.7 (native)"),
            "pocket_center": res.get("box_center"),
            "pocket_source": res.get("box_source", ""),
            "protein_pdb": protein_pdb,
            "target_name": target_name,
            "drug_name": drug_name,
            "note": ("Real AutoDock Vina 1.2.7 docking: receptor protonated (pdbfixer), "
                     "ligand prepared with Meeko, docked into the "
                     f"{res.get('box_source', 'binding site')}. ΔG in kcal/mol."),
        }
