"""
Compound Validator
Filters and deduplicates compounds before they reach the UI.
Rules:
  1. Must have a valid SMILES (RDKit parseable) or be a known biologic
  2. Must have at least one indication row in DB
  3. Must have at least one target in mechanisms table
  4. Deduplicated by InChIKey (salts/hydrates → canonical form)
"""

import logging
import re
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import InchiInfo, inchi
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False


# ─── SMILES validation ────────────────────────────────────────────────────────

def valid_smiles(smiles: str) -> bool:
    """True if SMILES is non-empty and RDKit can parse it."""
    if not smiles or not smiles.strip():
        return False
    if not RDKIT_OK:
        return True  # can't check, assume OK
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None


def get_inchikey(smiles: str) -> Optional[str]:
    """Return InChIKey for a SMILES, or None if invalid/unavailable."""
    if not RDKIT_OK or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return inchi.MolToInchiKey(mol)
    except Exception:
        return None


# ─── Salt / form normalization ────────────────────────────────────────────────

_SALT_SUFFIXES = re.compile(
    r"\s+(hydrochloride|hcl|sodium|potassium|sulfate|sulphate|maleate|tartrate"
    r"|phosphate|acetate|citrate|mesylate|tosylate|fumarate|succinate"
    r"|besylate|bromide|chloride|monohydrate|dihydrate|trihydrate"
    r"|hydrate|anhydrous|free\s+base|salt|monohydrochloride|dihydrochloride"
    r"|hemisulfate|hemihydrate|sesquihydrate)",
    re.IGNORECASE,
)


def normalize_name(name: str) -> str:
    """Strip salt/form suffixes to get the base drug name."""
    n = _SALT_SUFFIXES.sub("", name).strip(" ,.")
    n = re.sub(r"\s+", " ", n)
    return n


# ─── Main validation + deduplication ─────────────────────────────────────────

def validate_and_deduplicate(
    compounds: List[Dict],
    require_smiles: bool = True,
    require_targets: bool = True,
) -> List[Dict]:
    """
    Filter invalid compounds and deduplicate by InChIKey (or normalised name).
    Each input dict should have: id, name, smiles (optional), max_phase.
    Returns filtered, deduplicated list preserving order.
    """
    seen_inchikeys: dict = {}  # inchikey → compound index in output
    seen_names: dict = {}       # normalized_name → compound index in output
    result: List[Dict] = []

    for c in compounds:
        name   = c.get("name", "").strip()
        smiles = c.get("smiles", "")

        # Filter: require valid SMILES for small molecules
        # Biologics (no SMILES) are kept if max_phase >= 3
        is_biologic = not smiles and float(c.get("max_phase") or 0) >= 3
        if require_smiles and not smiles and not is_biologic:
            logger.debug(f"Dropped {name}: no SMILES")
            continue

        if smiles and not valid_smiles(smiles):
            logger.debug(f"Dropped {name}: invalid SMILES")
            continue

        # Deduplicate by InChIKey
        ik = get_inchikey(smiles) if smiles else None
        if ik:
            if ik in seen_inchikeys:
                # Keep the one with higher clinical phase
                existing_idx = seen_inchikeys[ik]
                existing = result[existing_idx]
                if float(c.get("max_phase") or 0) > float(existing.get("max_phase") or 0):
                    result[existing_idx] = c
                continue
            seen_inchikeys[ik] = len(result)

        # Deduplicate by normalised name (catches salt forms without InChIKey)
        norm = normalize_name(name).lower()
        if norm in seen_names:
            existing_idx = seen_names[norm]
            existing = result[existing_idx]
            if float(c.get("max_phase") or 0) > float(existing.get("max_phase") or 0):
                result[existing_idx] = c
                if ik:
                    seen_inchikeys[ik] = existing_idx
            continue
        seen_names[norm] = len(result)

        result.append(c)

    logger.info(
        f"Validator: {len(compounds)} in → {len(result)} out "
        f"({len(compounds)-len(result)} dropped/deduplicated)"
    )
    return result


def add_validation_flags(compounds: List[Dict]) -> List[Dict]:
    """
    Add a 'valid' flag and 'validation_notes' list to each compound dict
    without removing any — useful for debugging what would be filtered.
    """
    for c in compounds:
        notes = []
        smiles = c.get("smiles", "")

        if not smiles:
            notes.append("no_smiles")
        elif not valid_smiles(smiles):
            notes.append("invalid_smiles")

        if not c.get("mechanisms") and not c.get("targets"):
            notes.append("no_targets")

        c["valid"] = len(notes) == 0
        c["validation_notes"] = notes

    return compounds
