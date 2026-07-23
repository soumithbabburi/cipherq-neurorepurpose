"""
Lightweight physics-based molecular docking engine.

When no external docking engine (DiffDock / AutoDock Vina binary) is available, this
module performs genuine pocket-based docking rather than returning the ligand unmoved:

  1. Pocket detection   — centre on the co-crystallised ligand (the true binding site)
                          if present, else the most buried cavity centroid.
  2. Pose sampling      — RDKit conformers placed in the pocket with random rigid-body
                          orientations, each refined by Monte-Carlo hill-climbing.
  3. Scoring            — an AutoDock-Vina-style empirical function: attractive Gaussians
                          + steric repulsion + hydrophobic contact + hydrogen bonding,
                          divided by a torsional-entropy penalty → kcal/mol.

The ligand is actually translated/rotated into the binding site, so the returned SDF
poses render inside the protein and the affinities reflect shape/chemistry complementarity.
This is an approximation — clearly labelled as such — not a replacement for Vina/DiffDock.
"""

import logging
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    from scipy.spatial.distance import cdist as _cdist
except Exception:
    _cdist = None

logger = logging.getLogger(__name__)

# van der Waals radii (Å) — AutoDock-style heavy-atom radii
_VDW = {"C": 1.9, "N": 1.8, "O": 1.7, "S": 2.0, "P": 2.1,
        "F": 1.5, "CL": 1.8, "BR": 2.0, "I": 2.2, "H": 1.0}

# Vina-like scoring weights
_W_GAUSS1, _W_GAUSS2, _W_REPULSION = -0.0356, -0.00516, 0.840
_W_HYDROPHOBIC, _W_HBOND          = -0.0351, -0.587
_W_ROT                            = 0.0585

# Residue/atom names that are NOT a real bound ligand (waters, ions, buffers, sugars)
_NON_LIGAND = {
    "HOH", "WAT", "DOD", "H2O", "SO4", "PO4", "CL", "NA", "K", "MG", "CA",
    "ZN", "MN", "FE", "GOL", "EDO", "PEG", "ACT", "DMS", "TRS", "FMT", "NO3",
    "IOD", "BR", "CU", "NI", "CO", "CD", "MES", "EPE", "BME", "NAG", "MAN", "BMA",
}


def _vdw(elem: str) -> float:
    return _VDW.get(elem.upper(), 1.8)


# ── PDB parsing ──────────────────────────────────────────────────────────────

def _parse_protein(pdb_text: str):
    """Return (coords Nx3, radii N, is_donor N, is_acceptor N, is_hydrophobic N) for receptor heavy atoms."""
    coords, radii, don, acc, hyd = [], [], [], [], []
    het_by_res: Dict[str, list] = {}

    for ln in pdb_text.splitlines():
        rec = ln[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        elem = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if elem == "H":
            continue
        try:
            x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
        except ValueError:
            continue
        resname = ln[17:20].strip().upper()
        atomname = ln[12:16].strip().upper()

        if rec == "HETATM":
            # candidate co-crystal ligand atom (for pocket detection)
            if resname not in _NON_LIGAND:
                het_by_res.setdefault(f"{resname}_{ln[21:26]}", []).append((x, y, z))
            continue

        coords.append((x, y, z))
        radii.append(_vdw(elem))
        # Donor: N with H potential (backbone amide N, sidechain N). Acceptor: O (and some N).
        is_o = elem == "O"
        is_n = elem == "N"
        don.append(is_n)                       # N atoms can donate
        acc.append(is_o or atomname in ("ND1", "NE2"))  # O accepts; His N accepts
        hyd.append(elem in ("C", "S"))         # carbon/sulfur = hydrophobic

    if not coords:
        return None

    # Largest non-water HETATM group = bound ligand → pocket centre
    pocket_ligand = None
    if het_by_res:
        biggest = max(het_by_res.values(), key=len)
        if len(biggest) >= 6:                  # real drug-like ligand, not a stray ion
            pocket_ligand = np.array(biggest, dtype=float)

    return {
        "coords": np.array(coords, dtype=float),
        "radii":  np.array(radii, dtype=float),
        "don":    np.array(don, dtype=bool),
        "acc":    np.array(acc, dtype=bool),
        "hyd":    np.array(hyd, dtype=bool),
        "pocket_ligand": pocket_ligand,
    }


def _pocket_center(prot: dict) -> np.ndarray:
    """True binding site = bound-ligand centroid; else most buried protein atom neighbourhood."""
    if prot.get("pocket_ligand") is not None:
        return prot["pocket_ligand"].mean(axis=0)
    # Fallback: the atom with the most neighbours within 10 Å is the most buried → pocket-like
    coords = prot["coords"]
    if len(coords) > 1500:
        coords = coords[np.random.choice(len(coords), 1500, replace=False)]
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    buried = (d < 10.0).sum(axis=1)
    return coords[int(np.argmax(buried))]


# ── Ligand prep ──────────────────────────────────────────────────────────────

def _prep_ligand(smiles_or_sdf: str, n_conf: int = 6):
    """Return (rdkit_mol_with_Hs, heavy_idx, ligand atom-property arrays, n_rotatable)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors

    mol = None
    if "M  END" in smiles_or_sdf or smiles_or_sdf.lstrip().startswith(("\n", " ", "Mrv", "RDKit")):
        mol = Chem.MolFromMolBlock(smiles_or_sdf, removeHs=False)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles_or_sdf)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params))
    if not cids:
        if AllChem.EmbedMolecule(mol, params) != 0:
            return None
        cids = [0]
    for cid in cids:
        try:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid)
        except Exception:
            pass

    heavy_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() != "H"]
    props = _ligand_props(mol, heavy_idx)
    n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    return mol, cids, heavy_idx, props, n_rot


def _ligand_props(mol, heavy_idx) -> dict:
    """Per-heavy-atom donor/acceptor/hydrophobic flags and vdW radii."""
    don, acc, hyd, radii = [], [], [], []
    for i in heavy_idx:
        a = mol.GetAtomWithIdx(i)
        sym = a.GetSymbol().upper()
        is_n, is_o = sym == "N", sym == "O"
        n_h = a.GetTotalNumHs()
        don.append((is_n or is_o) and n_h > 0)          # donor if N/O bearing H
        acc.append(is_o or (is_n and a.GetTotalValence() < 4))
        hydrophobic = sym == "C" and all(
            nb.GetSymbol() in ("C", "H") for nb in a.GetNeighbors())
        hyd.append(hydrophobic)
        radii.append(_vdw(sym))
    return {
        "don": np.array(don, dtype=bool),
        "acc": np.array(acc, dtype=bool),
        "hyd": np.array(hyd, dtype=bool),
        "radii": np.array(radii, dtype=float),
    }


def score_pose(protein_pdb: str, pose_sdf: str) -> Optional[float]:
    """Empirical binding-energy estimate (kcal/mol) for an already-positioned pose
    (e.g. a DiffDock output). Returns None if it can't be scored."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        prot = _parse_protein(protein_pdb)
        if prot is None or len(prot["coords"]) < 8:
            return None
        mol = Chem.MolFromMolBlock(pose_sdf, removeHs=False)
        if mol is None:
            return None
        conf = mol.GetConformer()
        heavy_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() != "H"]
        coords = conf.GetPositions()[heavy_idx]
        lp = _ligand_props(mol, heavy_idx)
        n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
        center = coords.mean(axis=0)                      # pose sits in the pocket already
        Psub = _pocket_subset(prot, center)
        return round(float(_score(coords, lp, Psub, n_rot)), 2)
    except Exception as e:
        logger.debug(f"score_pose failed: {e}")
        return None


# ── Geometry helpers ─────────────────────────────────────────────────────────

def _rand_rotation(rng) -> np.ndarray:
    """Uniform random rotation matrix from a random quaternion."""
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([
        [1 - 2*(y*y+z*z), 2*(x*y-z*w),     2*(x*z+y*w)],
        [2*(x*y+z*w),     1 - 2*(x*x+z*z), 2*(y*z-x*w)],
        [2*(x*z-y*w),     2*(y*z+x*w),     1 - 2*(x*x+y*y)],
    ])


def _small_rotation(rng, sigma=0.35) -> np.ndarray:
    ax = rng.normal(size=3)
    ax /= (np.linalg.norm(ax) + 1e-9)
    ang = rng.normal() * sigma
    c, s = np.cos(ang), np.sin(ang)
    x, y, z = ax
    return np.array([
        [c+x*x*(1-c),   x*y*(1-c)-z*s, x*z*(1-c)+y*s],
        [y*x*(1-c)+z*s, c+y*y*(1-c),   y*z*(1-c)-x*s],
        [z*x*(1-c)-y*s, z*y*(1-c)+x*s, c+z*z*(1-c)],
    ])


# ── Scoring ──────────────────────────────────────────────────────────────────

def _score(lig_heavy: np.ndarray, lp: dict, P: dict, n_rot: int) -> float:
    """AutoDock-Vina-style empirical score (kcal/mol; lower = better).

    Only atom pairs within 8 Å contribute (every term is ~0 beyond that), so we
    evaluate the expensive Gaussians/clip ONLY on those pairs instead of the full
    ligand×receptor matrix. Identical result to the dense form, far less work.
    """
    pc = P["coords"]
    if _cdist is not None:
        d = _cdist(lig_heavy, pc)                    # (m, n) — C-optimised
    else:
        diff = lig_heavy[:, None, :] - pc[None, :, :]
        d = np.sqrt((diff * diff).sum(axis=2))
    wi, wj = np.nonzero(d < 8.0)
    if wi.size == 0:
        return 0.0
    t = d[wi, wj] - (lp["radii"][wi] + P["radii"][wj])   # surface distance, 1-D

    g1 = np.exp(-(t / 0.5) ** 2)
    g2 = np.exp(-((t - 3.0) / 2.0) ** 2)
    # Soft-core repulsion: a rigid conformer can't relax torsional strain the way a
    # flexible-ligand engine does, so cap each clash contribution (atomic softness).
    rep = np.where(t < 0, np.minimum(t * t, 3.0), 0.0)

    hyd_pair = lp["hyd"][wi] & P["hyd"][wj]
    hydrophobic = np.clip((1.5 - t) / 1.0, 0.0, 1.0) * hyd_pair

    hb_pair = (lp["don"][wi] & P["acc"][wj]) | (lp["acc"][wi] & P["don"][wj])
    hbond = np.clip((0.0 - t) / 0.7, 0.0, 1.0) * hb_pair

    inter = (_W_GAUSS1 * g1.sum()
             + _W_GAUSS2 * g2.sum()
             + _W_REPULSION * rep.sum()
             + _W_HYDROPHOBIC * hydrophobic.sum()
             + _W_HBOND * hbond.sum())
    return inter / (1.0 + _W_ROT * n_rot)


def _pocket_subset(P: dict, center: np.ndarray, radius: float = 14.0) -> dict:
    """Restrict receptor atoms to the pocket neighbourhood for speed."""
    d = np.linalg.norm(P["coords"] - center, axis=1)
    keep = d < radius
    if keep.sum() < 10:
        keep = d < (radius + 10.0)
    return {k: P[k][keep] for k in ("coords", "radii", "don", "acc", "hyd")}


# ── Main entry ───────────────────────────────────────────────────────────────

def dock_ligand(protein_pdb: str, ligand_smiles_or_sdf: str, drug_name: str = "ligand",
                n_poses: int = 5, seed: int = 7, center=None) -> Dict:
    """Dock a ligand into a protein pocket. Returns poses (SDF) + affinities (kcal/mol).

    If `center` (x,y,z) is given, docking targets that site — pass a named pocket's
    centre so the docked pose lands in the same pocket the UI labels/highlights.
    """
    try:
        from rdkit import Chem
    except Exception:
        return {"success": False, "error": "RDKit required for docking", "poses": []}

    prot = _parse_protein(protein_pdb)
    if prot is None or len(prot["coords"]) < 8:
        return {"success": False, "error": "No usable receptor atoms in PDB", "poses": []}

    prepared = _prep_ligand(ligand_smiles_or_sdf)
    if prepared is None:
        return {"success": False, "error": "Could not build 3-D ligand", "poses": []}
    mol, cids, heavy_idx, lp, n_rot = prepared

    if center is not None:
        center = np.asarray(center, dtype=float)
        pocket_source = "named binding pocket"
    else:
        center = _pocket_center(prot)
        pocket_source = "co-crystallised ligand" if prot.get("pocket_ligand") is not None else "largest buried cavity"
    Psub = _pocket_subset(prot, center)

    rng = np.random.default_rng(seed)
    n_trials = 16            # random placements per conformer
    n_mc     = 70            # simulated-annealing steps per placement
    box      = 3.0           # Å translational search half-width
    T0, T1   = 2.0, 0.10     # annealing temperature schedule

    candidates: List[Tuple[float, int, np.ndarray]] = []
    for cid in cids:
        conf = mol.GetConformer(cid)
        full = conf.GetPositions()                     # (N,3) incl. H
        heavy0 = full[heavy_idx]
        centroid = heavy0.mean(axis=0)
        base = full - centroid                         # centre the whole molecule

        for ti in range(n_trials):
            R = _rand_rotation(rng)
            # Seed several trials exactly on the pocket centre, others jittered
            trans = center + (np.zeros(3) if ti < 4 else rng.uniform(-box, box, size=3))
            cur = base @ R.T + trans                   # (N,3) placed
            cur_score = _score(cur[heavy_idx], lp, Psub, n_rot)
            best_cur, best_cur_score = cur, cur_score

            # Simulated-annealing rigid refinement (Metropolis acceptance)
            for step in range(n_mc):
                T = T0 * (T1 / T0) ** (step / max(1, n_mc - 1))
                cen = cur[heavy_idx].mean(axis=0)
                dR = _small_rotation(rng, sigma=0.45)
                dt = rng.normal(scale=0.9, size=3)
                cand = (cur - cen) @ dR.T + cen + dt
                s = _score(cand[heavy_idx], lp, Psub, n_rot)
                if s < cur_score or rng.random() < np.exp(-(s - cur_score) / T):
                    cur, cur_score = cand, s
                    if s < best_cur_score:
                        best_cur, best_cur_score = cand, s
            candidates.append((best_cur_score, cid, best_cur))

    if not candidates:
        return {"success": False, "error": "Docking produced no poses", "poses": []}

    candidates.sort(key=lambda c: c[0])

    # Diversify: keep poses that differ by > 2 Å heavy-atom RMSD
    chosen: List[Tuple[float, int, np.ndarray]] = []
    for sc, cid, coords in candidates:
        if len(chosen) >= n_poses:
            break
        h = coords[heavy_idx]
        if all(np.sqrt(((h - c2[heavy_idx]) ** 2).sum(axis=1).mean()) > 2.0 for _, _, c2 in chosen):
            chosen.append((sc, cid, coords))
    if not chosen:
        chosen = candidates[:n_poses]

    poses, affs, confs = [], [], []
    for sc, cid, coords in chosen:
        pose_mol = Chem.Mol(mol)
        conf = pose_mol.GetConformer(cid)
        for i in range(pose_mol.GetNumAtoms()):
            conf.SetAtomPosition(i, tuple(float(v) for v in coords[i]))
        try:
            pose_mol = Chem.RemoveHs(pose_mol)
        except Exception:
            pass
        poses.append(Chem.MolToMolBlock(pose_mol, confId=cid))
        aff = round(float(sc), 2)
        affs.append(aff)
        # Map affinity → 0–1 confidence (≈ -4 kcal/mol → low, ≈ -12 → high)
        confs.append(round(float(1.0 / (1.0 + np.exp((aff + 7.0) / 1.5))), 3))

    logger.info(f"dock_engine: {len(poses)} poses for {drug_name}; best {affs[0]} kcal/mol (pocket: {pocket_source})")
    return {
        "success": True,
        "poses": poses,
        "binding_affinities": affs,
        "confidence_scores": confs,
        "docking_method": "RepurposeIQ Dock (Vina-style empirical scoring)",
        "pocket_center": [round(float(v), 2) for v in center],
        "pocket_source": pocket_source,
        "note": ("Physics-based pocket docking (Vina-style empirical scoring): "
                 "attractive/repulsive vdW, hydrophobic contact, and H-bond terms. "
                 "An approximation for triage — confirm with AutoDock Vina/DiffDock."),
    }
