"""
Quantum Chemistry Engine — GFN2-xTB (semi-empirical QM)
======================================================

This is REAL quantum chemistry, not a force field. It runs Grimme's GFN2-xTB
semi-empirical tight-binding method (an approximate solution of the electronic
Schrödinger equation) on a drug-sized molecule in a few seconds, and returns
electronic-structure descriptors that medicinal chemists actually use:

  • HOMO / LUMO energies and the HOMO–LUMO gap (eV)   → reactivity, stability
  • Dipole moment (Debye), molecular polarizability    → permeability / solubility
  • QM atomic partial charges                          → electrostatics (better
                                                          than Gasteiger for docking)
  • Solvation free energy (ALPB implicit water, kcal)  → solubility / logP signal
  • Conceptual-DFT reactivity indices derived from the frontier orbitals
    (chemical hardness η, electronegativity χ, electrophilicity ω)

The xTB binary lives in an isolated micromamba env under .qc/ (created out of band),
so it does not touch the app's pip venv. RDKit (in the venv) builds the 3-D
geometry; xTB does the quantum calculation. Results are cached to disk (the
calculation is deterministic for a fixed seed).

If xTB is not installed, every entry point degrades cleanly to
{"available": False} — callers must handle that.
"""
import logging
import os
import re
import subprocess
import tempfile
from pathlib import Path
from shutil import which
from typing import Dict, List, Optional

from services.diskcache import DiskCache

logger = logging.getLogger(__name__)

_HARTREE_TO_KCAL = 627.5094740631
_ROOT = Path(__file__).parent.parent
_QC_DIR = _ROOT / ".qc"
_QC_ENV = os.getenv("QC_ENV", "qc")
_METHOD_LABEL = "GFN2-xTB (semi-empirical quantum mechanics)"

_cache = DiskCache("qc_props", ttl=None)  # deterministic → never expires
_micromamba_path: Optional[str] = None
_resolved = False


# ── Engine discovery ──────────────────────────────────────────────────────────
def _find_micromamba() -> Optional[str]:
    global _micromamba_path, _resolved
    if _resolved:
        return _micromamba_path
    _resolved = True
    env = os.getenv("QC_MICROMAMBA")
    candidates = [env] if env else []
    candidates.append(str(_QC_DIR / "Library" / "bin" / "micromamba.exe"))
    candidates.append(str(_QC_DIR / "Library" / "bin" / "micromamba"))
    for c in candidates:
        if c and Path(c).exists():
            _micromamba_path = c
            return c
    _micromamba_path = which("micromamba")
    return _micromamba_path


def available() -> bool:
    """True if the xTB engine can be invoked."""
    mm = _find_micromamba()
    if not mm:
        return False
    # The env dir existing is a good cheap proxy; a full --version call is avoided
    # on the hot path. Both win-native and unix layouts are checked.
    envdir = _QC_DIR / "envs" / _QC_ENV
    return envdir.exists()


def method_label() -> str:
    return _METHOD_LABEL


# ── xTB invocation ─────────────────────────────────────────────────────────────
def _run_xtb(xyz_text: str, charge: int, solvent: str, optimize: bool,
             timeout: int = 180) -> Dict:
    mm = _find_micromamba()
    if not mm:
        return {"available": False, "success": False, "error": "micromamba/xTB not found"}

    workdir = tempfile.mkdtemp(prefix="qc_")
    try:
        xyz_path = os.path.join(workdir, "mol.xyz")
        with open(xyz_path, "w", encoding="utf-8") as f:
            f.write(xyz_text)

        args = [mm, "-r", str(_QC_DIR), "run", "-n", _QC_ENV,
                "xtb", "mol.xyz", "--gfn", "2", "--chrg", str(charge),
                "--alpb", solvent, "--parallel", "2"]
        if optimize:
            args += ["--opt", "loose"]

        env = dict(os.environ)
        env["OMP_NUM_THREADS"] = "2"
        env["MKL_NUM_THREADS"] = "2"

        proc = subprocess.run(args, cwd=workdir, env=env, timeout=timeout,
                              capture_output=True, encoding="utf-8", errors="replace")
        out = proc.stdout or ""

        # xTB may error-terminate on a failed geometry opt yet still produce a valid
        # single-point block. Accept the result if the frontier orbitals were printed.
        parsed = _parse_xtb(out, workdir)
        if not parsed.get("homo_ev"):
            # Retry once as a pure single-point on the input geometry.
            if optimize:
                return _run_xtb(xyz_text, charge, solvent, optimize=False, timeout=timeout)
            return {"available": True, "success": False,
                    "error": "xTB produced no orbital data",
                    "stderr": (proc.stderr or "")[-400:]}
        parsed.update(available=True, success=True, method=_METHOD_LABEL,
                      solvent=solvent, charge=charge, optimized=optimize)
        return parsed
    except subprocess.TimeoutExpired:
        return {"available": True, "success": False, "error": "xTB timed out"}
    except Exception as e:
        logger.warning(f"xTB run failed: {e}")
        return {"available": True, "success": False, "error": str(e)}
    finally:
        try:
            import shutil
            shutil.rmtree(workdir, ignore_errors=True)
        except Exception:
            pass


def _parse_xtb(out: str, workdir: str) -> Dict:
    """Extract electronic-structure properties from xTB stdout + output files."""
    res: Dict = {}

    def _grab(pattern: str, group: int = 1) -> Optional[float]:
        m = re.search(pattern, out)
        if m:
            try:
                return float(m.group(group))
            except (ValueError, IndexError):
                return None
        return None

    # Frontier orbitals: lines end with "<energy_eV> (HOMO)" / "(LUMO)"
    for line in out.splitlines():
        if "(HOMO)" in line:
            toks = line.split()
            i = toks.index("(HOMO)")
            try:
                res["homo_ev"] = float(toks[i - 1])
            except (ValueError, IndexError):
                pass
        elif "(LUMO)" in line:
            toks = line.split()
            i = toks.index("(LUMO)")
            try:
                res["lumo_ev"] = float(toks[i - 1])
            except (ValueError, IndexError):
                pass

    gap = _grab(r"HOMO-LUMO GAP\s+([-\d.]+)\s*eV")
    if gap is not None:
        res["gap_ev"] = round(gap, 4)

    etot = _grab(r"::\s*total energy\s+([-\d.]+)\s*Eh")
    if etot is not None:
        res["total_energy_eh"] = round(etot, 6)

    gsolv = _grab(r"::\s*->\s*Gsolv\s+([-\d.]+)\s*Eh")
    if gsolv is not None:
        res["gsolv_kcal_mol"] = round(gsolv * _HARTREE_TO_KCAL, 3)

    # Dipole: the "full:" line inside the molecular-dipole block; last col = |μ| in D
    lines = out.splitlines()
    for i, ln in enumerate(lines):
        if "molecular dipole" in ln:
            for j in range(i, min(i + 5, len(lines))):
                if lines[j].strip().startswith("full:"):
                    toks = lines[j].split()
                    try:
                        res["dipole_debye"] = round(float(toks[-1]), 3)
                    except (ValueError, IndexError):
                        pass
            break

    alpha = _grab(r"Mol\.\s*α\(0\)\s*/au\s*:\s*([-\d.]+)")
    if alpha is None:
        alpha = _grab(r"Mol\.\s*alpha\(0\)\s*/au\s*:\s*([-\d.]+)")
    if alpha is not None:
        res["polarizability_au"] = round(alpha, 3)

    # Atomic partial charges (one per atom, heavy + H, in input order)
    charges_file = os.path.join(workdir, "charges")
    if os.path.exists(charges_file):
        try:
            vals = [float(x) for x in open(charges_file).read().split()]
            res["atomic_charges"] = [round(v, 4) for v in vals]
        except Exception:
            pass

    # Conceptual-DFT reactivity descriptors from the frontier orbitals
    homo, lumo = res.get("homo_ev"), res.get("lumo_ev")
    if homo is not None and lumo is not None:
        eta = (lumo - homo) / 2.0           # chemical hardness
        chi = -(homo + lumo) / 2.0          # electronegativity (Mulliken)
        res["hardness_ev"] = round(eta, 4)
        res["electronegativity_ev"] = round(chi, 4)
        res["chemical_potential_ev"] = round(-chi, 4)
        if eta > 1e-6:
            res["electrophilicity_ev"] = round((chi * chi) / (2.0 * eta), 4)

    return res


# ── Public API ─────────────────────────────────────────────────────────────────
def compute_properties(smiles: str, name: str = "", solvent: str = "water",
                       optimize: bool = False, use_cache: bool = True) -> Dict:
    """Compute GFN2-xTB electronic properties for a molecule from SMILES.

    Returns a dict with available/success flags plus the parsed descriptors.
    """
    if not available():
        return {"available": False, "success": False,
                "error": "Quantum engine (xTB) not installed", "method": _METHOD_LABEL}

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return {"available": True, "success": False, "error": "RDKit required for 3-D geometry"}

    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return {"available": True, "success": False, "error": "Invalid SMILES"}

    cache_key = f"{Chem.MolToSmiles(mol)}|{solvent}|opt={optimize}"
    if use_cache:
        cached = _cache.get(cache_key)
        if cached is not None:
            return cached

    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol_h, params) != 0:
        if AllChem.EmbedMolecule(mol_h, AllChem.ETDG()) != 0:
            return {"available": True, "success": False, "error": "3-D embedding failed"}
    try:
        AllChem.MMFFOptimizeMolecule(mol_h)
    except Exception:
        pass

    charge = Chem.GetFormalCharge(mol_h)
    conf = mol_h.GetConformer()
    lines = [str(mol_h.GetNumAtoms()), name or "molecule"]
    heavy_elems: List[str] = []
    for atom in mol_h.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        sym = atom.GetSymbol()
        lines.append(f"{sym} {p.x:.6f} {p.y:.6f} {p.z:.6f}")
        if sym != "H":
            heavy_elems.append(sym)
    xyz_text = "\n".join(lines) + "\n"

    result = _run_xtb(xyz_text, charge, solvent, optimize)
    if result.get("success"):
        result["n_atoms"] = mol_h.GetNumAtoms()
        result["formula"] = Chem.rdMolDescriptors.CalcMolFormula(mol)
        if use_cache:
            _cache.set(cache_key, result)
    return result


def qm_charges(smiles: str, solvent: str = "water") -> Optional[List[float]]:
    """Heavy-atom QM partial charges (RDKit heavy-atom order) for docking electrostatics,
    or None if unavailable. Hydrogens are dropped so the list aligns with dock_engine's
    heavy-atom arrays."""
    res = compute_properties(smiles, solvent=solvent)
    if not res.get("success") or "atomic_charges" not in res:
        return None
    try:
        from rdkit import Chem
        mol_h = Chem.AddHs(Chem.MolFromSmiles(smiles))
        charges = res["atomic_charges"]
        heavy = [charges[a.GetIdx()] for a in mol_h.GetAtoms() if a.GetSymbol() != "H"]
        return heavy
    except Exception:
        return None
