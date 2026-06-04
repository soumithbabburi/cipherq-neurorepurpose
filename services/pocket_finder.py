"""
Pocket Finder — detect, rank and NAME binding pockets on a protein structure.

Mirrors the NVIDIA-style docking output: instead of one anonymous binding site,
this surfaces the top-N cavities, scores their druggability, lines them with
their residues, and gives each a biological name where the structure tells us one
(co-crystallised ligand, PDB SITE records). Pure numpy + scipy — no fpocket binary.

Public API:
  find_pockets(pdb_text, max_pockets=5)  -> List[pocket dict]
  assign_pose_pocket(pose_sdf, pockets)  -> pocket id ("P1") or None

A pocket dict:
  { "id":"P1", "label":"Pocket 1", "name":"ATP-binding site",
    "center":[x,y,z], "volume": 320, "druggability": 0.78, "band":"High",
    "residues":["ASP228","THR231",...], "ligand":"ATP"|None, "n_points": 41 }
"""
import logging
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Waters / ions / buffers / cryo / sugars — never a real druggable ligand
_NON_LIGAND = {
    "HOH", "WAT", "DOD", "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN", "MN",
    "FE", "GOL", "EDO", "PEG", "PG4", "1PE", "ACT", "DMS", "TRS", "FMT", "NO3",
    "IOD", "BR", "CU", "NI", "CO", "CD", "MES", "EPE", "BME", "NAG", "MAN", "BMA",
    "FUC", "GAL", "BOG", "PGE", "MPD", "IMD", "AZI", "CIT", "TLA", "SCN",
}

# Hydrophobic residues — fraction lining a pocket is a strong druggability signal
_HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "PHE", "MET", "TRP", "PRO", "CYS"}

# Common cofactor / ligand HET codes → friendly biological site name
_LIGAND_SITE_NAME = {
    "ATP": "ATP-binding site", "ADP": "ADP-binding site", "AMP": "AMP-binding site",
    "ANP": "ATP-analog (AMP-PNP) site", "ACP": "ATP-analog (AMP-PCP) site",
    "GTP": "GTP-binding site", "GDP": "GDP-binding site", "GNP": "GTP-analog site",
    "NAD": "NAD(H) cofactor site", "NAI": "NADH cofactor site", "NAP": "NADP(H) site",
    "NDP": "NADPH site", "FAD": "FAD cofactor site", "FMN": "FMN cofactor site",
    "HEM": "Heme site", "HEC": "Heme C site", "SAM": "SAM methyl-donor site",
    "SAH": "SAH site", "COA": "Coenzyme-A site", "PLP": "PLP (B6) cofactor site",
    "TPP": "Thiamine-PP site", "BTN": "Biotin site", "RET": "Retinoid site",
}


# ── PDB parsing ──────────────────────────────────────────────────────────────
def _parse(pdb_text: str) -> Optional[dict]:
    """Parse protein heavy atoms, HET ligand groups, and SITE annotations."""
    coords: List[Tuple[float, float, float]] = []
    res_keys: List[str] = []        # "CHAIN:RESSEQ" per protein atom
    res_names: Dict[str, str] = {}  # res_key -> resname
    het_groups: Dict[str, dict] = {}
    site_residues: Dict[str, set] = {}     # site_id -> {res_key}
    site_desc: Dict[str, str] = {}         # site_id -> description
    _last_site_id = [None]

    for ln in pdb_text.splitlines():
        rec = ln[:6].strip()

        if rec == "REMARK" and ln[7:10].strip() == "800":
            body = ln[11:].strip()
            m = re.match(r"SITE_IDENTIFIER:\s*(\S+)", body)
            if m:
                _last_site_id[0] = m.group(1)
                continue
            m = re.match(r"SITE_DESCRIPTION:\s*(.+)", body)
            if m and _last_site_id[0]:
                site_desc[_last_site_id[0]] = m.group(1).strip()
            continue

        if rec == "SITE":
            sid = ln[11:14].strip()
            if not sid:
                continue
            bucket = site_residues.setdefault(sid, set())
            # up to 4 residues per record, fixed columns
            for off in (18, 29, 40, 51):
                rn = ln[off:off + 3].strip()
                ch = ln[off + 4:off + 5].strip()
                sq = ln[off + 5:off + 9].strip()
                if rn and sq:
                    bucket.add(f"{ch}:{sq}")
            continue

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
        chain = ln[21:22].strip()
        resseq = ln[22:26].strip()

        if rec == "HETATM":
            if resname not in _NON_LIGAND:
                key = f"{resname}_{chain}{resseq}"
                g = het_groups.setdefault(key, {"resname": resname, "chain": chain,
                                                "resseq": resseq, "coords": []})
                g["coords"].append((x, y, z))
            continue

        rk = f"{chain}:{resseq}"
        coords.append((x, y, z))
        res_keys.append(rk)
        res_names.setdefault(rk, resname)

    if len(coords) < 8:
        return None
    return {
        "coords": np.asarray(coords, dtype=float),
        "res_keys": res_keys,
        "res_names": res_names,
        "het_groups": het_groups,
        "site_residues": site_residues,
        "site_desc": site_desc,
    }


# ── Cavity detection ─────────────────────────────────────────────────────────
def _cavity_points(atoms: np.ndarray, tree, spacing: float) -> np.ndarray:
    """Grid points that sit in an enclosed cavity (empty, near surface, buried)."""
    lo = atoms.min(axis=0) - 4.0
    hi = atoms.max(axis=0) + 4.0
    # Cap total grid points for speed — coarsen spacing on very large structures
    n_est = np.prod((hi - lo) / spacing)
    while n_est > 260_000 and spacing < 2.4:
        spacing *= 1.25
        n_est = np.prod((hi - lo) / spacing)

    gx = np.arange(lo[0], hi[0], spacing)
    gy = np.arange(lo[1], hi[1], spacing)
    gz = np.arange(lo[2], hi[2], spacing)
    grid = np.array(np.meshgrid(gx, gy, gz, indexing="ij")).reshape(3, -1).T
    if not len(grid):
        return np.empty((0, 3))

    # 1. empty space near the surface (not inside protein, not far in solvent)
    d_near, _ = tree.query(grid, k=1)
    shell = grid[(d_near >= 2.6) & (d_near <= 4.6)]
    if not len(shell):
        return np.empty((0, 3))

    # 2. enclosure test — octant coverage by protein atoms within 9 Å
    keep = np.zeros(len(shell), dtype=bool)
    neigh = tree.query_ball_point(shell, r=9.0)
    for i, idxs in enumerate(neigh):
        if len(idxs) < 8:
            continue
        v = atoms[idxs] - shell[i]
        octs = ((v[:, 0] > 0).astype(int) * 4
                + (v[:, 1] > 0).astype(int) * 2
                + (v[:, 2] > 0).astype(int))
        if len(np.unique(octs)) >= 6:      # surrounded on most sides = buried cavity
            keep[i] = True
    return shell[keep]


def _cluster(points: np.ndarray, link: float) -> List[np.ndarray]:
    """Connected-components clustering of cavity points (single-linkage by distance)."""
    from scipy.spatial import cKDTree
    if not len(points):
        return []
    tree = cKDTree(points)
    pairs = tree.query_ball_tree(tree, r=link)
    seen = np.zeros(len(points), dtype=bool)
    clusters = []
    for start in range(len(points)):
        if seen[start]:
            continue
        stack, comp = [start], []
        seen[start] = True
        while stack:
            j = stack.pop()
            comp.append(j)
            for nb in pairs[j]:
                if not seen[nb]:
                    seen[nb] = True
                    stack.append(nb)
        clusters.append(points[comp])
    return clusters


def _druggability(volume: float, hydrophobic_frac: float, n_residues: int) -> float:
    """0–1 heuristic in the spirit of fpocket's drug score: size + lipophilicity + closure."""
    # volume sweet-spot ~150–600 Å³; saturate beyond
    vscore = max(0.0, min(1.0, (volume - 80.0) / 520.0))
    hscore = max(0.0, min(1.0, hydrophobic_frac / 0.55))
    rscore = max(0.0, min(1.0, (n_residues - 6) / 18.0))
    return round(0.45 * vscore + 0.35 * hscore + 0.20 * rscore, 2)


def _band(score: float) -> str:
    return "High" if score >= 0.6 else "Moderate" if score >= 0.35 else "Low"


def _clean_site_desc(desc: str) -> Optional[str]:
    """'BINDING SITE FOR RESIDUE ATP A 501' -> 'ATP-binding site'."""
    m = re.search(r"BINDING SITE FOR RESIDUE\s+([A-Z0-9]{1,3})", desc, re.I)
    if m:
        het = m.group(1).upper()
        return _LIGAND_SITE_NAME.get(het) or f"{het} binding site"
    m = re.search(r"BINDING SITE FOR (?:CHAIN|MONOMER)\s+(.+)", desc, re.I)
    if m:
        return f"Interface site ({m.group(1).strip()[:24]})"
    return desc.strip().capitalize()[:40] if desc else None


# ── Public API ───────────────────────────────────────────────────────────────
def find_pockets(pdb_text: str, max_pockets: int = 5) -> List[dict]:
    """Detect, rank and name up to `max_pockets` binding pockets."""
    try:
        from scipy.spatial import cKDTree
    except Exception:
        logger.warning("scipy unavailable — pocket detection skipped")
        return []
    if not pdb_text or len(pdb_text) < 200:
        return []
    P = _parse(pdb_text)
    if P is None:
        return []

    full_tree = cKDTree(P["coords"])
    spacing = 1.4

    def _lining(center: np.ndarray, radius: float = 5.0):
        """Residues + hydrophobic fraction lining a sphere at `center`."""
        idx = full_tree.query_ball_point(center, r=radius)
        res_keys = {P["res_keys"][i] for i in idx}
        ordered = sorted(res_keys, key=lambda k: (k.split(":")[0], int(k.split(":")[1] or 0)))
        names = [(rk, P["res_names"].get(rk, "")) for rk in ordered]
        hyd = (sum(1 for _, n in names if n in _HYDROPHOBIC_RES) / len(names)) if names else 0.0
        residues = [f"{n}{rk.split(':')[1]}" for rk, n in names if n][:12]
        return res_keys, residues, hyd

    site_res = P["site_residues"]
    site_desc = P["site_desc"]

    def _site_name(res_keys: set) -> Optional[str]:
        if not site_res:
            return None
        best_ov, best_sid = 0, None
        for sid, rks in site_res.items():
            ov = len(rks & res_keys)
            if ov > best_ov:
                best_ov, best_sid = ov, sid
        if best_sid and best_ov >= 2:
            return _clean_site_desc(site_desc.get(best_sid, "")) or f"Annotated site {best_sid}"
        return None

    pockets: List[dict] = []
    seed_centers: List[np.ndarray] = []

    # ── 1. Ligand-seeded pockets (highest confidence — a bound drug marks a site) ──
    for g in P["het_groups"].values():
        if len(g["coords"]) < 5:
            continue
        lig_atoms = np.asarray(g["coords"], dtype=float)
        center = lig_atoms.mean(axis=0)
        res_keys, residues, hyd = _lining(center, radius=5.5)
        if len(residues) < 4:
            continue
        # volume ≈ enclosed cavity around the ligand (ligand atom count is a fair proxy)
        volume = int(round(len(lig_atoms) * 22 + 60))
        drug = round(min(1.0, max(_druggability(volume, hyd, len(res_keys)), 0.72) + 0.05), 2)
        name = _LIGAND_SITE_NAME.get(g["resname"]) or _site_name(res_keys) or f"{g['resname']} binding site"
        # Collapse the same site duplicated across chains (A/B copies share residue names)
        sig = frozenset(residues)
        dup = next((q for q in pockets if q["ligand"] == g["resname"]
                    and len(sig & q["_sig"]) >= 0.6 * min(len(sig), len(q["_sig"]) or 1)), None)
        if dup:
            dup["copies"] = dup.get("copies", 1) + 1
            seed_centers.append(center)
            continue
        pockets.append({
            "center": [round(float(v), 2) for v in center], "_np": center, "_sig": sig,
            "volume": volume, "druggability": drug, "band": _band(drug),
            "residues": residues, "ligand": g["resname"], "bound": True,
            "name": name, "n_points": len(lig_atoms), "copies": 1,
        })
        seed_centers.append(center)

    # ── 2. Geometric (apo / allosteric) cavities for everything else ──────────────
    atoms = P["coords"]
    if len(atoms) > 24000:
        atoms = atoms[np.random.default_rng(0).choice(len(atoms), 24000, replace=False)]
    tree = cKDTree(atoms)
    clusters = _cluster(_cavity_points(atoms, tree, spacing), link=spacing * 1.1)
    geo = []
    for c in clusters:
        if len(c) < 6:
            continue
        center = c.mean(axis=0)
        volume = len(c) * (spacing ** 3)
        if volume < 80 or volume > 2200:        # discrete pocket, not a solvent channel
            continue
        if any(np.linalg.norm(center - s) < 7.0 for s in seed_centers):
            continue                            # already covered by a ligand seed
        res_keys, residues, hyd = _lining(center, radius=5.0)
        if len(residues) < 5:
            continue
        drug = _druggability(volume, hyd, len(res_keys))
        geo.append({
            "center": [round(float(v), 2) for v in center], "_np": center,
            "volume": int(round(volume)), "druggability": drug, "band": _band(drug),
            "residues": residues, "ligand": None, "bound": False, "copies": 1,
            "name": _site_name(res_keys), "n_points": int(len(c)),
        })
    geo.sort(key=lambda p: (p["druggability"], p["volume"]), reverse=True)

    # Dedup near-coincident geometric pockets (multi-chain copies)
    for p in geo:
        if any(np.linalg.norm(p["_np"] - q["_np"]) < 6.0 for q in pockets):
            continue
        pockets.append(p)

    pockets = pockets[:max_pockets]
    if not pockets:
        return []

    # Label P1..N (bound/ligand pockets already lead the list)
    for i, p in enumerate(pockets, start=1):
        p["id"] = f"P{i}"
        p["label"] = f"Pocket {i}"
        if not p.get("name"):
            p["name"] = p["label"]
        p.pop("_np", None)
        p.pop("_sig", None)

    logger.info("pocket_finder: %d pockets — %s", len(pockets),
                ", ".join(f"{p['label']}={p['name']}" for p in pockets))
    return pockets


def assign_pose_pocket(pose_sdf: str, pockets: List[dict], max_dist: float = 12.0) -> Optional[str]:
    """Return the id of the pocket nearest a docked pose's centroid (or None)."""
    if not pockets or not pose_sdf:
        return None
    try:
        coords = []
        in_atoms = False
        for ln in pose_sdf.splitlines():
            # V2000 counts line is line 4; atom block follows. Robust float-triple scan:
            parts = ln.split()
            if len(parts) >= 4:
                try:
                    x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                    # 4th token is the element symbol in an atom line
                    if parts[3].isalpha():
                        coords.append((x, y, z)); in_atoms = True
                        continue
                except ValueError:
                    pass
            if in_atoms and not coords:
                break
        if not coords:
            return None
        cen = np.mean(coords, axis=0)
        best_id, best_d = None, max_dist
        for p in pockets:
            d = float(np.linalg.norm(np.asarray(p["center"]) - cen))
            if d < best_d:
                best_d, best_id = d, p["id"]
        return best_id
    except Exception:
        return None
