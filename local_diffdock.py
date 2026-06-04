"""
Local DiffDock Integration
Runs DiffDock from a local installation (https://github.com/gcorso/DiffDock).

Setup:
  git clone https://github.com/gcorso/DiffDock
  cd DiffDock && pip install -e .
  Set env var: DIFFDOCK_PATH=/path/to/DiffDock

Usage:
  runner = LocalDiffDock()
  if runner.available:
      result = runner.run(protein_pdb, ligand_sdf, drug_name)
"""

import logging
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


class LocalDiffDock:
    """Runs DiffDock locally via subprocess."""

    def __init__(self):
        self.diffdock_path = self._find_diffdock()
        self.python_exec = self._find_python()
        self.available = self._check_available()
        if self.available:
            logger.info(f"LocalDiffDock ready at {self.diffdock_path}")
        else:
            logger.info("LocalDiffDock not found — set DIFFDOCK_PATH env var to enable")

    # ── Discovery ─────────────────────────────────────────────────────────────
    def _find_diffdock(self) -> str:
        """Find DiffDock installation from env var or common paths."""
        env = os.environ.get("DIFFDOCK_PATH", "")
        if env and Path(env, "inference.py").exists():
            return env

        candidates = [
            Path.home() / "DiffDock",
            Path.home() / "diffdock",
            Path("DiffDock"),
            Path("diffdock"),
            Path.cwd().parent / "DiffDock",
        ]
        for p in candidates:
            if (p / "inference.py").exists():
                return str(p)
        return ""

    def _find_python(self) -> str:
        for name in ("python", "python3"):
            if shutil.which(name):
                return name
        return "python"

    def _check_available(self) -> bool:
        if not self.diffdock_path:
            return False
        inf = Path(self.diffdock_path) / "inference.py"
        return inf.exists()

    # ── Main entry ────────────────────────────────────────────────────────────
    def run(
        self,
        protein_pdb: str,
        ligand_sdf: str,
        drug_name: str = "ligand",
        num_poses: int = 10,
        timeout: int = 600,
    ) -> Dict:
        """
        Run DiffDock inference and return poses + scores.

        Returns dict with keys:
            success (bool), poses (list[str]), confidence_scores (list[float]),
            binding_affinities (list[float]), docking_method (str), error (str)
        """
        if not self.available:
            return {"success": False, "error": "DiffDock not installed", "poses": []}

        with tempfile.TemporaryDirectory(prefix="diffdock_") as tmpdir:
            tmp = Path(tmpdir)
            prot_file = tmp / "protein.pdb"
            lig_file  = tmp / "ligand.sdf"
            out_dir   = tmp / "output"
            out_dir.mkdir()

            # Clean PDB: keep only ATOM records of the first/largest chain.
            # Raw RCSB PDBs contain water (HOH), HETATM co-crystals, ANISOU records
            # and sometimes NaN coordinates — all of which cause DiffDock's SVD to fail.
            cleaned = self._clean_pdb(protein_pdb)
            prot_file.write_text(cleaned)
            lig_file.write_text(ligand_sdf)

            cmd = [
                self.python_exec,
                str(Path(self.diffdock_path) / "inference.py"),
                "--protein_path",   str(prot_file),
                "--ligand_description", str(lig_file),
                "--out_dir",        str(out_dir),
                "--inference_steps","20",
                "--samples_per_complex", str(num_poses),
                "--batch_size",     "10",
                "--no_final_step_noise",
            ]

            logger.info(f"Running DiffDock: {' '.join(cmd[:4])}… (protein: {len(cleaned)} bytes)")
            try:
                proc = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                    cwd=self.diffdock_path,
                )
            except subprocess.TimeoutExpired:
                return {"success": False, "error": f"DiffDock timed out after {timeout}s", "poses": []}
            except Exception as e:
                return {"success": False, "error": str(e), "poses": []}

            if proc.returncode != 0:
                err = (proc.stderr or proc.stdout or "")[-600:]
                return {"success": False, "error": f"DiffDock exit {proc.returncode}: {err}", "poses": []}

            result = self._parse_output(out_dir, drug_name)
            # Attach cleaned PDB so caller can use it for the 3-D viewer
            if result.get("success"):
                result["protein_pdb"] = cleaned
                # DiffDock outputs a *confidence* (log-odds of a good pose), NOT a
                # binding free energy. Keep confidence as the pose-quality metric and
                # estimate a real ΔG per pose with the empirical scorer.
                try:
                    from services.dock_engine import score_pose
                    scored = [score_pose(cleaned, p) for p in result["poses"]]
                    if any(s is not None for s in scored):
                        conf = result.get("confidence_scores", [])
                        affs = [
                            s if s is not None else (round(-conf[i] * 1.36, 2) if i < len(conf) else 0.0)
                            for i, s in enumerate(scored)
                        ]
                        # Re-rank poses by empirical binding energy (best ΔG first),
                        # keeping each pose's DiffDock confidence paired with it.
                        order = sorted(range(len(affs)), key=lambda i: affs[i])
                        result["poses"]              = [result["poses"][i] for i in order]
                        result["binding_affinities"] = [affs[i] for i in order]
                        result["confidence_scores"]  = [conf[i] if i < len(conf) else 0.0 for i in order]
                        result["affinity_method"]    = "Empirical (Vina-style) scoring of DiffDock pose"
                except Exception as e:
                    logger.debug(f"DiffDock pose rescore failed: {e}")
            return result

    # ── PDB preprocessing ────────────────────────────────────────────────────
    @staticmethod
    def _clean_pdb(pdb_text: str) -> str:
        """
        Strip raw RCSB PDB down to protein ATOM records only.
        Removes: HETATM (waters, ions, co-crystals), ANISOU, REMARK, SEQRES etc.
        Keeps only the largest chain (by residue count) to avoid multi-chain SVD issues.
        """
        lines = pdb_text.splitlines()

        # Collect ATOM records, skip water residues (HOH / WAT / DOD)
        WATER = {"HOH", "WAT", "DOD", "H2O"}
        atom_lines: List[str] = [
            ln for ln in lines
            if ln.startswith("ATOM") and ln[17:20].strip() not in WATER
        ]

        if not atom_lines:
            # Fallback: accept HETATM too (some structures use them for protein)
            atom_lines = [ln for ln in lines if ln.startswith(("ATOM", "HETATM"))]

        # Count residues per chain (column 21 = chain ID in standard PDB)
        chain_counts: dict = {}
        for ln in atom_lines:
            chain = ln[21] if len(ln) > 21 else "A"
            resnum = ln[22:26].strip() if len(ln) > 26 else "0"
            chain_counts.setdefault(chain, set()).add(resnum)

        # Pick largest chain
        if chain_counts:
            best_chain = max(chain_counts, key=lambda c: len(chain_counts[c]))
            atom_lines = [ln for ln in atom_lines if len(ln) > 21 and ln[21] == best_chain]

        # Remove duplicate atoms (keep first occurrence per atom serial)
        seen: set = set()
        unique: List[str] = []
        for ln in atom_lines:
            serial = ln[6:11].strip()
            if serial not in seen:
                seen.add(serial)
                unique.append(ln)

        if not unique:
            return pdb_text  # give up: return original if nothing survived

        return "\n".join(unique) + "\nEND\n"

    # ── Output parsing ────────────────────────────────────────────────────────
    def _parse_output(self, out_dir: Path, drug_name: str) -> Dict:
        """Parse DiffDock output directory into standard result dict."""
        poses: List[str]     = []
        confs: List[float]   = []
        affs:  List[float]   = []

        # DiffDock writes: rank1_confidence-1.23.sdf, rank2_confidence0.45.sdf …
        # Sort by NUMERIC rank (lexicographic sort would give rank1, rank10, rank2 …).
        def _rank_num(p):
            m = re.search(r"rank(\d+)", p.name)
            return int(m.group(1)) if m else 999
        sdf_files = sorted(out_dir.rglob("rank*_confidence*.sdf"), key=_rank_num)
        if not sdf_files:
            sdf_files = sorted(out_dir.rglob("*.sdf"))

        for sdf_path in sdf_files:
            try:
                poses.append(sdf_path.read_text())
            except Exception:
                continue

            # Parse confidence from filename: rank1_confidence-1.23.sdf
            conf = self._parse_confidence(sdf_path.name)
            confs.append(conf)
            # Convert confidence (log-odds-like) to kcal/mol approximation
            affs.append(round(-conf * 1.36, 2))

        if not poses:
            return {"success": False, "error": "DiffDock produced no poses", "poses": []}

        logger.info(f"DiffDock: {len(poses)} poses for {drug_name}")
        return {
            "success": True,
            "poses": poses,
            "confidence_scores": confs,
            "binding_affinities": affs,
            "docking_method": "Local DiffDock",
        }

    @staticmethod
    def _parse_confidence(filename: str) -> float:
        """Extract confidence value from DiffDock output filename."""
        import re
        # Use non-greedy decimal so trailing '.' before '.sdf' is not captured
        m = re.search(r"confidence(-?\d+(?:\.\d+)?)", filename)
        return float(m.group(1)) if m else 0.0

    # ── Installation helper ───────────────────────────────────────────────────
    @staticmethod
    def install_instructions() -> str:
        return (
            "To enable local DiffDock:\n"
            "  1. git clone https://github.com/gcorso/DiffDock\n"
            "  2. cd DiffDock && pip install -e .\n"
            "  3. Download model weights (see DiffDock README)\n"
            "  4. Set env var: DIFFDOCK_PATH=/absolute/path/to/DiffDock\n"
            "  5. Restart the app"
        )
