"""
Drug modality classifier — route each asset to the right optimization workflow.
═══════════════════════════════════════════════════════════════════════════════
SMILES/atom-level optimization only makes sense for SMALL MOLECULES. An antibody
(mepolizumab, adalimumab…) must be optimized in sequence/developability space, not
by editing atoms. This classifier detects the modality so the platform can branch:

    small molecule → SMILES workshop (scaffold edits, ADMET, docking)
    antibody       → sequence-based biologics workflow (variants + developability)
    peptide/protein/oligonucleotide → biologic handling (no SMILES optimization)

Signals (most reliable first): explicit ChEMBL molecule_type; INN stem (-mab = mAb,
-cept = Fc-fusion, -tide = peptide, -sen/-mersen = antisense oligo); presence/size of
a parseable SMILES. Fail-soft → 'small_molecule' only when a real small-molecule SMILES
is present, else 'unknown'.
"""
from __future__ import annotations

import re
from typing import Dict, Optional

# INN suffix stems (WHO nomenclature) → modality.
_MAB = re.compile(r"(mab|mab\b)$", re.I)            # monoclonal antibody
_ANTIBODY_STEMS = ("mab",)
_FUSION_STEMS = ("cept",)                            # etanercept, aflibercept (Fc-fusion)
_PEPTIDE_STEMS = ("tide", "relin", "actide", "pressin", "eptide")
_OLIGO_STEMS = ("sen", "mersen", "rsen", "virsen")   # antisense/siRNA (-sen)
_PROTEIN_HINTS = ("interferon", "insulin", "erythropoietin", "epoetin", "filgrastim",
                  "somatropin", "factor viii", "factor ix", "albumin", "enzyme")


def classify(name: str = "", smiles: str = "", molecule_type: str = "",
             compound: Optional[Dict] = None) -> Dict:
    """Return {modality, optimization, is_small_molecule, confidence, basis}.

    modality ∈ small_molecule | antibody | fusion_protein | peptide | protein |
              oligonucleotide | unknown
    optimization ∈ 'smiles' | 'sequence'
    """
    if compound:
        name = name or compound.get("name", "")
        smiles = smiles or compound.get("smiles", "")
        molecule_type = molecule_type or compound.get("molecule_type", "")
    n = (name or "").strip().lower()
    mt = (molecule_type or "").strip().lower()
    smi = (smiles or "").strip()

    def out(modality, opt, conf, basis):
        return {"modality": modality, "optimization": opt,
                "is_small_molecule": modality == "small_molecule",
                "confidence": conf, "basis": basis,
                "label": modality.replace("_", " ")}

    # 1. Explicit ChEMBL molecule_type (most authoritative)
    if mt:
        if "antibody" in mt:
            return out("antibody", "sequence", "high", "ChEMBL molecule_type = Antibody")
        if "protein" in mt:
            return out("protein", "sequence", "high", "ChEMBL molecule_type = Protein")
        if "oligonucleotide" in mt or "oligo" in mt:
            return out("oligonucleotide", "sequence", "high", "ChEMBL molecule_type")
        if "enzyme" in mt or "gene" in mt or "cell" in mt:
            return out("protein", "sequence", "high", f"ChEMBL molecule_type = {mt}")
        if "small molecule" in mt:
            return out("small_molecule", "smiles", "high", "ChEMBL molecule_type = Small molecule")

    # 2. INN stem on the name (works even with no molecule_type)
    if n.endswith(_ANTIBODY_STEMS):
        return out("antibody", "sequence", "high", "INN stem '-mab' (monoclonal antibody)")
    if n.endswith(_FUSION_STEMS):
        return out("fusion_protein", "sequence", "medium", "INN stem '-cept' (Fc-fusion protein)")
    if n.endswith(_OLIGO_STEMS):
        return out("oligonucleotide", "sequence", "medium", "INN stem (antisense/oligo)")
    if any(h in n for h in _PROTEIN_HINTS):
        return out("protein", "sequence", "medium", "recognized therapeutic-protein name")
    if n.endswith(_PEPTIDE_STEMS):
        return out("peptide", "sequence", "medium", "INN stem (peptide)")

    # 3. SMILES-based: a real, non-trivial small-molecule structure
    if smi:
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                mw = Descriptors.MolWt(mol)
                # very large SMILES (>1500 Da) is usually a peptide/biologic drawn out
                if mw > 1500 and mol.GetNumAtoms() > 100:
                    return out("peptide", "sequence", "low",
                               f"large structure (MW {mw:.0f}) — likely a peptide/biologic")
                return out("small_molecule", "smiles", "high", "parseable small-molecule SMILES")
        except Exception:
            pass
        return out("small_molecule", "smiles", "low", "SMILES present (RDKit unavailable)")

    # 4. No SMILES and no biologic signal → unknown (do NOT assume small molecule)
    return out("unknown", "sequence", "low",
               "no SMILES and no small-molecule evidence — modality undetermined")


if __name__ == "__main__":
    import json
    for nm, smi, mt in [
        ("mepolizumab", "", ""),
        ("adalimumab", "", "Antibody"),
        ("imatinib", "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1", "Small molecule"),
        ("etanercept", "", ""),
        ("semaglutide", "", ""),
        ("nusinersen", "", ""),
    ]:
        print(f"{nm:14} -> {json.dumps(classify(nm, smi, mt))}")
