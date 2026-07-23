"""
AutoDock Vina docking worker  —  runs under the `drugdisc-agent` conda env.
═══════════════════════════════════════════════════════════════════════════════
This is executed as a SUBPROCESS by services/vina_engine.py (which runs in the
app's own venv). All the heavy cheminformatics deps (RDKit, Meeko, pdbfixer,
OpenBabel) live only in the drugdisc-agent env, so the app never has to import
them — it just shells out to this worker and reads back JSON.

Full native-Windows chain (no Vina Python API, which has no win-64 build):
  1. Ligand : SMILES -> RDKit 3-D (ETKDGv3 + MMFF) -> Meeko PDBQT
  2. Receptor: PDB text -> strip heterogens/waters -> pdbfixer (missing atoms +
               hydrogens at pH 7.0) -> OpenBabel rigid PDBQT
  3. Box     : caller-supplied center, else the co-crystal ligand centroid (the
               REAL binding site), else the receptor centroid (blind, flagged)
  4. Dock    : the official AutoDock Vina Windows CLI binary (vina.exe)
  5. Poses   : parse `REMARK VINA RESULT` affinities, split models, convert each
               to SDF (OpenBabel) so the existing 3Dmol viewer can render them

Contract: reads a JSON job on argv[1] (a file path), writes a JSON result to the
--out path. Never raises to the shell without also emitting a JSON error.
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

# Residue names that are NOT a druggable co-crystal ligand (solvent/ions/buffers).
_NON_LIGAND = {
    "HOH", "WAT", "DOD", "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN", "MN",
    "FE", "GOL", "EDO", "PEG", "ACT", "DMS", "TRS", "FMT", "NO3", "IOD", "BR",
    "CU", "NI", "CO", "CD", "MES", "EPE", "BME", "NAG", "MAN", "BMA", "PLP",
}


def _log(msg: str):
    print(msg, file=sys.stderr, flush=True)


# ── Receptor preparation ─────────────────────────────────────────────────────

def _cocrystal_ligand_centroid(pdb_text: str):
    """Centroid of the largest drug-like HETATM group = the true binding site."""
    groups = {}
    for ln in pdb_text.splitlines():
        if ln[:6].strip() != "HETATM":
            continue
        res = ln[17:20].strip().upper()
        if res in _NON_LIGAND:
            continue
        try:
            x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
        except ValueError:
            continue
        groups.setdefault(f"{res}_{ln[21:26]}", []).append((x, y, z))
    if not groups:
        return None
    biggest = max(groups.values(), key=len)
    if len(biggest) < 6:                       # too small to be a real ligand
        return None
    n = len(biggest)
    return [round(sum(c[i] for c in biggest) / n, 3) for i in range(3)]


def _prepare_receptor(pdb_text: str, work: Path) -> Path:
    """Clean + protonate the receptor and write a rigid PDBQT for Vina."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    raw = work / "receptor_raw.pdb"
    raw.write_text(pdb_text, encoding="utf-8")

    fixer = PDBFixer(filename=str(raw))
    fixer.removeHeterogens(keepWater=False)     # drop ligand/ions/waters -> apo receptor
    fixer.findMissingResidues()
    # Only rebuild missing atoms in existing residues, not whole missing loops
    # (which pdbfixer would model unreliably and can create phantom geometry).
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    # Do NOT rebuild missing heavy atoms. OpenMM's addMissingAtoms() segfaults
    # (native access violation, rc=0xC0000005) on some structures — e.g. the ABL1 and
    # KIT kinase domains — which kills the whole worker and silently drops docking to a
    # weak physics fallback (the "binding affinities not working" bug). Missing side-
    # chain atoms are almost always surface, not pocket, so skipping the rebuild costs
    # little for docking; the hydrogens that actually matter for scoring are still added.
    fixer.missingAtoms = {}
    fixer.missingTerminals = {}
    fixer.addMissingHydrogens(7.0)

    fixed = work / "receptor_fixed.pdb"
    with open(fixed, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh)

    # PDB -> rigid PDBQT via OpenBabel (Gasteiger charges, receptor mode = -xr)
    from openbabel import pybel
    mol = next(pybel.readfile("pdb", str(fixed)))
    rec_pdbqt = work / "receptor.pdbqt"
    mol.write("pdbqt", str(rec_pdbqt), overwrite=True, opt={"r": True})
    # OpenBabel sometimes emits ROOT/branch records for a receptor; strip anything
    # that isn't an atom so Vina reads it as a pure rigid receptor.
    lines = [ln for ln in rec_pdbqt.read_text().splitlines()
             if ln[:4] == "ATOM" or ln[:6] == "HETATM"]
    if len(lines) < 20:
        raise RuntimeError("receptor PDBQT has too few atoms after prep")
    rec_pdbqt.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return rec_pdbqt


# ── Ligand preparation ───────────────────────────────────────────────────────

def _prepare_ligand(smiles: str, work: Path) -> Path:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"RDKit could not parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        # retry with random coords for awkward molecules
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise RuntimeError("RDKit failed to embed a 3-D conformer")
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass

    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups:
        raise RuntimeError("Meeko produced no MoleculeSetup")
    pdbqt_str, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError(f"Meeko PDBQT write failed: {err}")
    lig = work / "ligand.pdbqt"
    lig.write_text(pdbqt_str, encoding="utf-8")
    return lig


# ── Vina + pose conversion ───────────────────────────────────────────────────

def _run_vina(vina_exe: str, rec: Path, lig: Path, center, size,
              exhaustiveness: int, n_poses: int, work: Path):
    out = work / "docked.pdbqt"
    cmd = [
        vina_exe,
        "--receptor", str(rec), "--ligand", str(lig),
        "--center_x", str(center[0]), "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness), "--num_modes", str(n_poses),
        "--seed", "42", "--out", str(out),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if proc.returncode != 0 or not out.exists():
        raise RuntimeError(f"vina.exe failed (rc={proc.returncode}): "
                           f"{(proc.stderr or proc.stdout)[-500:]}")
    return out


def _split_and_convert(out_pdbqt: Path, work: Path):
    """Split the multi-model docked PDBQT into per-pose SDF blocks + affinities."""
    from openbabel import pybel

    text = out_pdbqt.read_text(encoding="utf-8", errors="ignore")
    models, cur, affs = [], [], []
    for ln in text.splitlines():
        if ln.startswith("MODEL"):
            cur = []
        elif ln.startswith("REMARK VINA RESULT"):
            try:
                affs.append(float(ln.split()[3]))
            except (IndexError, ValueError):
                affs.append(None)
            cur.append(ln)
        elif ln.startswith("ENDMDL"):
            models.append("\n".join(cur))
        else:
            cur.append(ln)
    if not models:                              # single model, no MODEL records
        models = [text]

    sdf_poses = []
    for i, block in enumerate(models):
        pth = work / f"pose_{i}.pdbqt"
        pth.write_text(block if block.strip().startswith(("MODEL", "REMARK", "ATOM", "HETATM"))
                       else block, encoding="utf-8")
        try:
            m = next(pybel.readfile("pdbqt", str(pth)))
            sdf_poses.append(m.write("sdf"))
        except Exception as e:
            _log(f"pose {i} SDF conversion failed: {e}")
    return sdf_poses, affs


def main():
    job = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    out_path = Path(job["out_json"])
    result = {"success": False}
    try:
        work = Path(tempfile.mkdtemp(prefix="vina_"))
        pdb_text = Path(job["protein_pdb_file"]).read_text(encoding="utf-8", errors="ignore")

        # Binding site: explicit center wins; else the co-crystal ligand (real
        # pocket); else the receptor centroid (blind docking, flagged honestly).
        center = job.get("center")
        box_source = "caller-supplied center"
        if not center:
            center = _cocrystal_ligand_centroid(pdb_text)
            box_source = "co-crystal ligand (orthosteric site)"
        if not center:
            box_source = "receptor centroid (blind — no ligand/site given)"

        rec = _prepare_receptor(pdb_text, work)
        if not center:                          # fall back to receptor centroid
            from openbabel import pybel
            m = next(pybel.readfile("pdbqt", str(rec)))
            cs = [a.coords for a in m.atoms]
            center = [round(sum(c[i] for c in cs) / len(cs), 3) for i in range(3)]

        lig = _prepare_ligand(job["smiles"], work)
        size = job.get("size", [22.0, 22.0, 22.0])
        out_pdbqt = _run_vina(job["vina_exe"], rec, lig, center, size,
                              job.get("exhaustiveness", 8), job.get("n_poses", 9), work)
        poses, affs = _split_and_convert(out_pdbqt, work)
        if not poses:
            raise RuntimeError("docking produced no parseable poses")

        result = {
            "success": True,
            "poses": poses,
            "binding_affinities": [a for a in affs if a is not None][:len(poses)],
            "box_center": [round(float(c), 3) for c in center],
            "box_size": size,
            "box_source": box_source,
            "docking_method": "AutoDock Vina 1.2.7 (native)",
        }
    except Exception as e:
        import traceback
        result = {"success": False, "error": str(e), "trace": traceback.format_exc()[-1200:]}
    out_path.write_text(json.dumps(result), encoding="utf-8")


if __name__ == "__main__":
    main()
