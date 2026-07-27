"""
Docking Service — priority chain:
  1. Local DiffDock (downloadable, free)
  2. NVIDIA DiffDock API (requires NVIDIA_API_KEY)
  3. AutoDock Vina (RDKit-based fallback)
"""
import logging
import os
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


class DockingService:
    def __init__(self, disease_name: str = ""):
        self.disease_name = disease_name

        # Try local DiffDock first
        try:
            from local_diffdock import LocalDiffDock
            self._local = LocalDiffDock()
        except Exception:
            self._local = None

        # NVIDIA DiffDock (API) — only attempt if local DiffDock is not available
        self._nvidia = None
        self._nvidia_ok = False
        if not (self._local and self._local.available):
            try:
                from nvidia_bionemo_integration import NVIDIABioNeMoClient
                self._nvidia = NVIDIABioNeMoClient()
                self._nvidia_ok = self._nvidia.api_available
            except Exception:
                pass

        # Boltz-2 binding-affinity engine (P4) — opt-in, GPU/weights required. Off
        # the default path; used only when method="boltz" AND the CLI is installed.
        self._boltz = None
        try:
            from services.boltz_engine import BoltzEngine
            _b = BoltzEngine()
            self._boltz = _b if _b.available else None
        except Exception:
            self._boltz = None

        # Real AutoDock Vina (native Windows binary + drugdisc-agent env). When
        # present this is the DEFAULT engine — it docks the ligand into the real
        # binding site with proper receptor prep and returns a true ΔG, replacing
        # the built-in Vina-STYLE approximation for the 'fast' path. Fail-soft.
        self._vina = None
        try:
            from services.vina_engine import VinaEngine
            _v = VinaEngine()
            self._vina = _v if _v.available else None
        except Exception:
            self._vina = None

        # Always available: the built-in physics-based pocket engine (RDKit + numpy)
        # is the guaranteed fallback, so docking always returns real spatial poses.
        self.available = True

        method = (
            "AutoDock Vina 1.2.7 (native)" if self._vina
            else "Local DiffDock" if (self._local and self._local.available)
            else "NVIDIA DiffDock" if self._nvidia_ok
            else "RepurposeIQ Dock (physics-based pocket engine)"
        )
        logger.info(f"DockingService ready — primary: {method}")

    # ── Public API ─────────────────────────────────────────────────────────────
    def perform_docking(
        self,
        drug_name: str,
        target_name: str,
        ligand_smiles: Optional[str] = None,
        protein_pdb: Optional[str] = None,
        method: str = "fast",
    ) -> Dict:
        smiles = ligand_smiles or self._get_drug_smiles(drug_name)
        if not smiles:
            return {"success": False, "error": f"No SMILES for {drug_name}", "poses": []}

        ligand_sdf = self._smiles_to_sdf(smiles, drug_name)
        if not ligand_sdf:
            return {"success": False, "error": "SMILES→SDF conversion failed (RDKit required)", "poses": []}

        # Use pre-fetched protein if available; otherwise fetch (slower)
        if not protein_pdb or len(protein_pdb) < 200:
            protein_pdb = self._get_protein_structure(target_name)

        # Tier 4 — failure honesty: if no real structure was found we fall back to a
        # generic placeholder peptide. Don't dress that up as a real result.
        no_structure = self._is_placeholder_structure(protein_pdb)

        # Detect & name pockets up front so the physics engine can dock INTO the
        # functional pocket (keeps the pose, its residues and the pocket label all
        # coincident in the viewer, like a DiffDock result).
        pockets = []
        if not no_structure:
            try:
                from services.pocket_finder import find_pockets
                pockets = find_pockets(protein_pdb, max_pockets=5)
            except Exception as e:
                logger.warning(f"pocket detection failed: {e}")

        # Tier 2 — pocket relevance: prefer a pocket lined by a co-crystallised ligand
        # (a real, functional binding site) over the merely-largest geometric cavity.
        dock_center = None
        if pockets:
            functional = next((p for p in pockets if p.get("bound")), pockets[0])
            dock_center = functional["center"]

        result = self._dock(protein_pdb, smiles, ligand_sdf, drug_name, target_name, method, dock_center)
        # Tag each pose with the pocket it landed in (reuse the pockets found above)
        self._attach_pockets(result, protein_pdb, pockets)
        # Tier 2 — validate poses (clash + pocket containment) on the receptor used
        self._validate_poses(result, protein_pdb)
        # Tier 3 — consensus: cross-score externally-generated poses (DiffDock/Vina)
        # with our empirical function so agreement can be judged.
        if not no_structure:
            self._add_consensus(result, protein_pdb)

        if no_structure and isinstance(result, dict):
            result["structure_warning"] = (
                f"No experimental or predicted structure was found for {target_name}. "
                "Showing a ligand-only view — no protein docking was performed.")
            result["structure_quality"] = "none"
            result["protein_pdb"] = ""   # ligand-only render rather than a fake blob
        elif isinstance(result, dict):
            result.setdefault("structure_quality", "real")
        return result

    @staticmethod
    def _is_placeholder_structure(pdb: str) -> bool:
        """True if the receptor is the generic fallback peptide (or effectively empty)."""
        if not pdb or len(pdb) < 200:
            return True
        if "GENERIC PROTEIN" in pdb:
            return True
        # A real receptor has many residues; the placeholder has a handful.
        return pdb.count("ATOM ") < 40

    def _add_consensus(self, result: Dict, protein_pdb: str) -> None:
        """Cross-score externally generated poses (DiffDock / Vina) with the empirical
        function for a consensus signal. Skipped for our own empirical engine, whose
        affinities already come from that function."""
        if not isinstance(result, dict) or not result.get("success"):
            return
        method = (result.get("docking_method") or "").lower()
        if "cipherq" in method or "empirical" in method:
            return  # already empirically scored
        poses = result.get("poses") or []
        if not poses:
            return
        try:
            from services.dock_engine import score_pose
            scores = [score_pose(protein_pdb, p) for p in poses]
            if any(s is not None for s in scores):
                result["empirical_affinities"] = scores
                result["consensus_note"] = (
                    "Poses re-scored with the empirical (Vina-style) function as an "
                    "independent consensus check on the deep-learning ranking.")
        except Exception as e:
            logger.debug(f"consensus scoring skipped: {e}")

    def _validate_poses(self, result: Dict, protein_pdb: str) -> None:
        """Annotate each pose with clash/containment validity; drop clearly bad poses
        when better ones exist. A no-op if there's nothing spatial to check."""
        if not isinstance(result, dict) or not result.get("success"):
            return
        poses = result.get("poses") or []
        if not poses or not protein_pdb:
            return
        try:
            import numpy as np
            from services.dock_engine import _parse_protein
            prot = _parse_protein(protein_pdb)
            if prot is None or len(prot["coords"]) < 8:
                return
            pcoords = prot["coords"]

            def _pose_atoms(sdf: str):
                pts, in_atoms = [], False
                for ln in sdf.splitlines():
                    parts = ln.split()
                    if len(parts) >= 4:
                        try:
                            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                            if parts[3].isalpha():
                                pts.append((x, y, z)); in_atoms = True
                                continue
                        except ValueError:
                            pass
                    if in_atoms and not pts:
                        break
                return np.asarray(pts, dtype=float) if pts else None

            flags = []
            for sdf in poses:
                lig = _pose_atoms(sdf)
                if lig is None or not len(lig):
                    flags.append({"valid": False, "reason": "no coords"}); continue
                # nearest receptor atom per ligand atom
                d = np.linalg.norm(lig[:, None, :] - pcoords[None, :, :], axis=2).min(axis=1)
                severe_clashes = int((d < 1.6).sum())          # hard steric overlap
                contacts = int((d < 4.5).sum())                # is it nestled in the protein?
                contained = contacts >= max(3, 0.3 * len(lig))
                valid = contained and severe_clashes <= 1
                flags.append({"valid": bool(valid),
                              "clashes": severe_clashes,
                              "contained": bool(contained)})
            result["pose_valid"] = flags
            # If some poses are valid and others aren't, keep only the valid ones
            # (preserving rank order) so the UI never shows a floating/clashing pose.
            valid_idx = [i for i, fl in enumerate(flags) if fl.get("valid")]
            if valid_idx and len(valid_idx) < len(poses):
                keep = set(valid_idx)
                for listkey in ("poses", "confidence_scores", "binding_affinities", "pose_pockets"):
                    seq = result.get(listkey)
                    if isinstance(seq, list) and len(seq) == len(poses):
                        result[listkey] = [v for i, v in enumerate(seq) if i in keep]
                result["pose_valid"] = [fl for i, fl in enumerate(flags) if i in keep]
        except Exception as e:
            logger.debug(f"pose validation skipped: {e}")

    def _dock(self, protein_pdb, smiles, ligand_sdf, drug_name, target_name, method="fast", dock_center=None) -> Dict:
        # Boltz-2 (opt-in): predicts a BINDING AFFINITY, not just a pose. Requires the
        # installed CLI (+GPU); falls back to the physics engine if unavailable or it errors.
        if method == "boltz":
            if self._boltz:
                result = self._boltz.run(protein_pdb, smiles, drug_name, target_name)
                if result.get("success"):
                    return result
                logger.warning(f"Boltz-2 failed: {result.get('error')} — falling back to physics engine")
            else:
                logger.info("Boltz-2 requested but not installed — using physics engine")
            return self._run_fallback(protein_pdb, smiles, ligand_sdf, drug_name, target_name, dock_center)

        # Default path ("fast" or explicit "vina"): prefer REAL AutoDock Vina when
        # the engine is available and we have a real receptor. It docks into the
        # detected pocket (dock_center) and returns a true ΔG. Falls back to the
        # physics engine if Vina is unavailable, the structure is a placeholder, or
        # the run fails — so behaviour degrades gracefully, never breaks.
        if method in ("fast", "vina"):
            if self._vina and not self._is_placeholder_structure(protein_pdb):
                vres = self._vina.run(protein_pdb, smiles, drug_name, target_name,
                                      dock_center=list(dock_center) if dock_center is not None else None)
                if vres.get("success"):
                    return vres
                logger.warning(f"Vina failed: {vres.get('error')} — falling back to physics engine")
            return self._run_fallback(protein_pdb, smiles, ligand_sdf, drug_name, target_name, dock_center)

        # Any other non-DiffDock method also uses the physics engine.
        if method != "diffdock":
            return self._run_fallback(protein_pdb, smiles, ligand_sdf, drug_name, target_name, dock_center)

        # 1. Local DiffDock (opt-in — slow but higher accuracy)
        if self._local and self._local.available:
            result = self._local.run(protein_pdb, ligand_sdf, drug_name)
            if result.get("success"):
                result.update(target_name=target_name, drug_name=drug_name)
                return result
            logger.warning(f"Local DiffDock failed: {result.get('error')} — trying next")

        # 2. NVIDIA DiffDock
        if self._nvidia_ok:
            result = self._nvidia.run_diffdock(
                protein_pdb=protein_pdb,
                ligand_sdf=ligand_sdf,
                ligand_name=drug_name,
                num_poses=20,
            )
            if result.get("success") and result.get("status") != "failed":
                result.update(
                    docking_method="NVIDIA DiffDock",
                    target_name=target_name,
                    drug_name=drug_name,
                )
                return result
            logger.warning(f"NVIDIA DiffDock failed: {result.get('error')} — falling back to physics engine")

        # Fall back to the built-in physics-based pocket docking engine.
        return self._run_fallback(protein_pdb, smiles, ligand_sdf, drug_name, target_name, dock_center)

    def _attach_pockets(self, result: Dict, protein_pdb: str, pockets=None) -> None:
        """Tag each pose with the pocket it landed in. Reuses pockets found earlier
        (or detects them now if not supplied)."""
        if not isinstance(result, dict) or not result.get("success"):
            return
        pdb = result.get("protein_pdb") or protein_pdb
        if not pdb or len(pdb) < 200:
            return
        try:
            from services.pocket_finder import find_pockets, assign_pose_pocket
            if pockets is None:
                pockets = find_pockets(pdb, max_pockets=5)
            if not pockets:
                return
            result["pockets"] = pockets
            poses = result.get("poses", []) or []
            result["pose_pockets"] = [assign_pose_pocket(p, pockets) for p in poses]
        except Exception as e:
            logger.warning(f"pocket detection failed: {e}")

    def is_available(self) -> bool:
        # The physics-based pocket engine is always available (RDKit + numpy),
        # so docking is never reduced to a non-spatial descriptor estimate.
        return True

    # ── Fallback chain ────────────────────────────────────────────────────────
    def _run_fallback(self, protein_pdb, smiles, ligand_sdf, drug_name, target_name, dock_center=None) -> Dict:
        # Built-in physics-based pocket docking (real pose generation + scoring).
        try:
            from services.dock_engine import dock_ligand
            res = dock_ligand(protein_pdb, smiles or ligand_sdf, drug_name, n_poses=5, center=dock_center)
            if res.get("success"):
                res.update(target_name=target_name, drug_name=drug_name)
                res["protein_pdb"] = protein_pdb   # same frame as the docked poses
                return res
            logger.warning(f"dock_engine: {res.get('error')}")
        except Exception as e:
            logger.error(f"dock_engine failed: {e}")

        # 3c. Last resort: descriptor estimate (no spatial pose)
        return self._rdkit_score(ligand_sdf, drug_name, target_name)

    def _rdkit_score(self, ligand_sdf, drug_name, target_name) -> Dict:
        """Generate estimated scores using RDKit descriptors when no docking engine is available."""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, AllChem
            mol = Chem.MolFromMolBlock(ligand_sdf)
            if mol is None:
                raise ValueError("Could not parse SDF")
            mw   = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            hbd  = Descriptors.NumHDonors(mol)
            hba  = Descriptors.NumHAcceptors(mol)
            # Property-based ΔG estimate — NOT a docked pose. Return a single,
            # honest number: never fabricate multiple "poses" (the same ligand
            # repeated) or invented confidence scores. This path only runs if the
            # physics engine itself failed, which is why there is no 3-D pose here.
            aff = round(-(0.8 + logp*0.3 - tpsa*0.01 + hbd*0.2 + hba*0.15), 2)
            # This descriptor regression is NOT a docking ΔG — no pose was computed.
            # Keep it OUT of `binding_affinities` so it can never be rendered as a
            # "Docking ΔG (kcal/mol)"; expose it only in a distinctly-named field. The
            # docking-ΔG slot stays empty (honest: nothing was docked on this fallback).
            return {
                "success": True,
                "poses": [ligand_sdf],            # ligand geometry only — not docked
                "binding_affinities": [],
                "descriptor_estimate": aff,       # property-based, no physical ΔG meaning
                "docking_method": "RDKit descriptor estimate (no 3-D pose)",
                "target_name": target_name,
                "drug_name": drug_name,
                "structure_quality": "none",
                "note": "No protein pose could be generated for this target; only a rough "
                        "property-based estimate is available (not a docking ΔG).",
            }
        except Exception as e:
            return {"success": False, "error": f"All docking methods failed: {e}", "poses": []}

    # ── Helpers ───────────────────────────────────────────────────────────────
    def _get_protein_structure(self, target_name: str) -> str:
        try:
            from real_pdb_fetcher import RealPDBFetcher
            fetcher = RealPDBFetcher()
            pdb = fetcher.fetch_protein_structure(target_name)
            if pdb and len(pdb) > 200:
                return pdb
        except Exception as e:
            logger.debug(f"PDB fetch failed for {target_name}: {e}")
        return self._generic_protein()

    def _get_drug_smiles(self, drug_name: str) -> Optional[str]:
        """Look up a drug's SMILES from the compounds table (canonical DB config).

        Matches on name (case-insensitive) or ChEMBL id. Returns None if not found
        so the caller can report an honest 'No SMILES' error rather than docking a
        wrong molecule.
        """
        if not drug_name:
            return None
        try:
            import psycopg2
            from config import db_params
            conn = psycopg2.connect(**db_params())
            cur = conn.cursor()
            cur.execute(
                "SELECT smiles FROM compounds "
                "WHERE (name ILIKE %s OR chembl_id = %s) AND smiles IS NOT NULL "
                "ORDER BY name = %s DESC LIMIT 1",
                (drug_name, drug_name, drug_name),
            )
            row = cur.fetchone()
            cur.close(); conn.close()
            return row[0] if row else None
        except Exception as e:
            logger.debug(f"SMILES lookup failed for {drug_name}: {e}")
            return None

    def _smiles_to_sdf(self, smiles: str, name: str = "molecule") -> Optional[str]:
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            mol = Chem.MolFromSmiles(smiles)
            if mol is None: return None
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            mol.SetProp("_Name", name)
            return Chem.MolToMolBlock(mol)
        except Exception as e:
            logger.error(f"SMILES→SDF failed: {e}")
            return None

    @staticmethod
    def _generic_protein() -> str:
        return """HEADER    GENERIC PROTEIN
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.267   2.395   0.000  1.00 20.00           O
ATOM      5  CB  ALA A   1       1.993  -0.742   1.228  1.00 20.00           C
ATOM      6  N   GLY A   2       3.331   1.520   0.000  1.00 20.00           N
ATOM      7  CA  GLY A   2       4.025   2.802   0.000  1.00 20.00           C
ATOM      8  C   GLY A   2       5.532   2.621   0.000  1.00 20.00           C
ATOM      9  O   GLY A   2       6.061   1.511   0.000  1.00 20.00           O
ATOM     10  N   LEU A   3       6.241   3.746   0.000  1.00 20.00           N
ATOM     11  CA  LEU A   3       7.694   3.746   0.000  1.00 20.00           C
ATOM     12  C   LEU A   3       8.246   5.166   0.000  1.00 20.00           C
ATOM     13  O   LEU A   3       7.504   6.141   0.000  1.00 20.00           O
ATOM     14  CB  LEU A   3       8.229   2.999   1.228  1.00 20.00           C
END
"""
