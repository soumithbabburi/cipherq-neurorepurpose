"""
Medicinal-chemistry advisor — WHAT to actually modify, and why.
================================================================
The old optimization panel gave generic, property-threshold advice ("PSA high →
reduce PSA for BBB"). That is shallow: it doesn't see the *real* liabilities a
med-chemist optimizes against — a metabolically labile ester, a mutagenic nitro-
aromatic, a phenol that will be glucuronidated to nothing, a Michael acceptor that
binds covalently, an aniline that bioactivates. This module reads the actual
structure and returns a PRIORITIZED liability analysis: each finding names the
offending substructure, WHY it's a problem, the SPECIFIC modification, and the
trade-off — then identifies the single limiting liability.

Signals combined:
  1. Structural alerts   — RDKit FilterCatalog (PAINS / Brenk / NIH)
  2. Metabolic soft spots & reactive groups — curated SMARTS with real fixes
  3. Property liabilities — CNS-MPO limiting term (services/cns_mpo) + area-aware
     developability out-of-range props (services/developability)

Everything is bounded and honest: no structure → "unavailable". RDKit optional.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
    from rdkit.Chem import FilterCatalog
    from rdkit.Chem.FilterCatalog import FilterCatalogParams
    _RDKIT = True
except Exception:
    _RDKIT = False

_SEV_RANK = {"high": 3, "medium": 2, "low": 1}

# ── Curated structure→liability rules ────────────────────────────────────────
# (name, detection SMARTS, severity, why, specific modification, trade-off, rxn).
# `rxn` is an RDKit reaction SMARTS that APPLIES the recommended fix in one click
# (empty when the fix is contextual and can't be auto-applied — advice only then).
# These encode real ADME-Tox knowledge, not physicochemical thresholds.
_LIABILITY_RULES = [
    ("Aromatic nitro group", "[c][$([NX3+](=O)[O-]),$([NX3](=O)=O)]", "high",
     "Nitroaromatics are a classic mutagenicity/genotoxicity alert — nitroreductases "
     "generate reactive nitroso/hydroxylamino metabolites that adduct DNA.",
     "Remove the nitro or replace it with a non-genotoxic isostere (nitrile, carboxamide, "
     "or a small EWG like -CF3 / -SO2Me) that keeps the electron-withdrawing effect.",
     "May shift electronics/potency at the binding site — re-check the SAR at that position.",
     "[c:1][N+](=O)[O-]>>[c:1]"),

    # Reactive Michael acceptor ONLY — the β-carbon (far from the carbonyl) must NOT bear
    # N/O/S. That exclusion drops stabilized push-pull enaminones / vinylogous amides
    # (e.g. the 3-(aminomethylene)oxindole core of nintedanib/sunitinib), which are
    # electronically deactivated and NOT electrophilic — a common false positive.
    ("Michael acceptor (α,β-unsaturated carbonyl)",
     "[CX3;!$([CX3][#7,#8,#16]):1]=[CX3][CX3]=[OX1]", "high",
     "An electrophilic Michael acceptor reacts covalently with thiols (glutathione, off-"
     "target cysteines) — idiosyncratic toxicity and rapid GSH clearance unless the "
     "covalent warhead is intended.",
     "Saturate the double bond, or add β-substitution / steric bulk to blunt reactivity; "
     "if covalency IS the design, tune the warhead's intrinsic reactivity instead.",
     "Saturation can lower potency if the planar enone is part of the pharmacophore.",
     "[CX3;!$([CX3][#7,#8,#16]):1]=[CX3:2][CX3:3]=[OX1:4]>>[C:1][C:2][C:3]=[O:4]"),

    ("Primary aromatic amine (aniline)", "[NX3;H2][c]", "high",
     "Anilines bioactivate (CYP/N-oxidation, acetylation) to hydroxylamine/nitroso species "
     "linked to mutagenicity and methemoglobinemia.",
     "Cap the amine (acylate to an amide, or N-alkylate), or block the para/ortho positions "
     "to stop the bioactivation pathway.",
     "Acylation adds MW/PSA and removes an H-bond donor — may change target engagement.",
     "[NX3;H2:1][c:2]>>[c:2][N:1]C(C)=O"),

    ("Alkyl/benzylic halide", "[CX4][Cl,Br,I]", "high",
     "A reactive alkyl/benzylic halide is a direct-acting alkylator (DNA/protein) — a "
     "genotoxicity liability and chemically unstable.",
     "Replace the leaving group with a stable substituent (–F if a small EWG is needed, "
     "–OMe, –Me), or move to a non-benzylic position.",
     "Loss of a synthetic handle; electronics change.",
     "[CX4:1][Cl,Br,I]>>[C:1]"),

    ("Epoxide / aziridine", "[OX2r3,NX3r3]1[#6][#6]1", "high",
     "Strained three-membered rings are potent electrophiles — alkylate nucleophiles "
     "(genotoxic, reactive-metabolite risk).",
     "Open the ring to a diol/amino-alcohol, or replace with a non-strained isostere.",
     "Changes 3-D shape and H-bonding substantially.",
     ""),

    ("Ester (metabolic soft spot)", "[CX3](=[OX1])[OX2H0][#6]", "medium",
     "Esters are cleaved rapidly by ubiquitous serum/tissue esterases → short half-life "
     "(unless a prodrug is intended).",
     "Replace the ester with a metabolically stable bioisostere — amide, 1,2,4-oxadiazole, "
     "or oxazole — to block esterase hydrolysis while preserving H-bond geometry.",
     "Amides raise PSA/HBD (can hurt permeability/BBB); oxadiazoles keep it flatter.",
     "[CX3:1](=[OX1:2])[OX2H0][#6]>>[C:1](=[O:2])N"),

    ("Phenol (phase-II conjugation site)", "[c][OX2H]", "medium",
     "Free phenols are rapidly glucuronidated/sulfated → high clearance and low oral "
     "exposure.",
     "Methylate the phenol, or use a bioisostere (F, or an intramolecular H-bond mask) if "
     "the –OH is not essential for target H-bonding.",
     "The phenol may be a key H-bond donor to the target — verify before masking.",
     "[c:1][OX2H1:2]>>[c:1][O:2]C"),

    ("Furan / thiophene (bioactivation-prone)", "[$(c1ccoc1),$(c1ccsc1)]", "medium",
     "Furans and thiophenes are epoxidized/S-oxidized by CYPs to reactive metabolites "
     "(hepatotoxicity signal).",
     "Swap for a more metabolically robust heteroaromatic — thiazole, oxazole, or a "
     "substituted phenyl — at that vector.",
     "Ring geometry/electronics differ; re-optimize the vector.",
     ""),

    ("Thiourea / thioamide", "[NX3][CX3]=[SX1]", "medium",
     "Thioureas/thioamides are metabolically activated (S-oxidation) and are frequent "
     "tox/PAINS flags.",
     "Replace C=S with C=O (urea/amide) or a cyanoguanidine isostere.",
     "Reduces H-bond strength/geometry — re-check potency.",
     "[NX3:1][CX3:2]=[SX1]>>[N:1][C:2]=O"),

    ("Hydrazine / hydrazide", "[NX3][NX3]", "medium",
     "Hydrazines are genotoxic and hepatotoxic (reactive diazonium/radical metabolites).",
     "Replace with an amide/carbamate linkage or a stable N-heterocycle.",
     "Backbone geometry change.",
     ""),

    ("Aldehyde", "[CX3H1](=O)[#6]", "medium",
     "Aldehydes are reactive electrophiles (Schiff-base formation with lysines) and are "
     "quickly oxidized/reduced in vivo.",
     "Oxidize to the acid/amide, reduce to the alcohol, or mask as a heterocycle if the "
     "carbonyl is not needed.",
     "Removes an electrophilic contact — potency-dependent.",
     "[CX3H1:1]=[OX1:2]>>[C:1](=[O:2])O"),
]


def _limiting_property_advice(mw, logp, tpsa, hbd, hba, rtb, arom) -> List[Dict]:
    """Property-level liabilities framed by the mechanism they hurt (not just 'BBB')."""
    out: List[Dict] = []
    if logp is not None and logp > 5:
        out.append({
            "type": "property", "severity": "high", "property": "logp",
            "title": f"High lipophilicity (cLogP {logp:.1f})",
            "detail": "cLogP > 5 drives poor aqueous solubility, high plasma-protein binding, "
                      "CYP metabolism and promiscuity/hERG risk — the dominant liability here.",
            "action": "Lower cLogP by ~1–2 units: introduce a polar handle (–OH, –NH2, morpholine) "
                      "at a solvent-exposed vector, or replace a lipophilic aryl with a polar heteroaryl.",
            "tradeoff": "Added polarity can cut permeability/BBB — place it away from the pharmacophore."})
    elif logp is not None and logp < 1:
        out.append({
            "type": "property", "severity": "medium", "property": "logp",
            "title": f"Low lipophilicity (cLogP {logp:.1f})",
            "detail": "cLogP < 1 usually means weak membrane permeability and low cell/CNS uptake.",
            "action": "Add measured lipophilicity — a methyl/ethyl, F, or replace a polar group with a "
                      "lipophilic bioisostere.",
            "tradeoff": "Watch solubility and off-target promiscuity as logP rises."})
    if tpsa is not None and tpsa > 140:
        out.append({
            "type": "property", "severity": "medium", "property": "tpsa",
            "title": f"High polar surface area (TPSA {tpsa:.0f} Å²)",
            "detail": "TPSA > 140 Å² predicts poor passive permeability and essentially no oral/CNS "
                      "absorption.",
            "action": "Trim H-bonding groups or mask a donor (N-methylation, intramolecular H-bond) to "
                      "bring TPSA under ~120 Å² (or ~90 for CNS).",
            "tradeoff": "Removing donors/acceptors can weaken target binding."})
    if hbd is not None and hbd > 5:
        out.append({
            "type": "property", "severity": "medium", "property": "hbd",
            "title": f"Too many H-bond donors ({hbd})",
            "detail": "HBD > 5 hurts passive permeability and desolvation penalty on binding.",
            "action": "Cap or bioisosterically replace donors — N-methylate an amide N–H, or convert "
                      "an –OH donor to an ether/F.",
            "tradeoff": "Each donor removed may be a target contact — verify the H-bond map."})
    if rtb is not None and rtb > 10:
        out.append({
            "type": "property", "severity": "low", "property": "rtb",
            "title": f"High flexibility ({rtb} rotatable bonds)",
            "detail": "Excessive rotatable bonds raise the conformational entropy penalty on binding "
                      "and lower oral bioavailability.",
            "action": "Rigidify — cyclize a flexible chain or add a ring/macrocyclic constraint that "
                      "locks the bioactive conformation.",
            "tradeoff": "Wrong constraint geometry can abolish binding — needs the bound pose."})
    if arom is not None and arom > 3:
        out.append({
            "type": "property", "severity": "low", "property": "arom_rings",
            "title": f"Many aromatic rings ({arom})",
            "detail": "High aromatic-ring count correlates with poor solubility, higher hERG and "
                      "developability attrition (Ritchie/Macdonald).",
            "action": "Replace one aromatic ring with a saturated/partially-saturated bioisostere "
                      "(e.g. phenyl → bicyclo[1.1.1]pentane or a saturated heterocycle).",
            "tradeoff": "sp3 replacements change vector geometry and may cost potency."})
    return out


def _structural_alerts(mol) -> List[Dict]:
    out: List[Dict] = []
    try:
        params = FilterCatalogParams()
        for cat in ("PAINS", "BRENK", "NIH"):
            params.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, cat))
        catalog = FilterCatalog.FilterCatalog(params)
        seen = set()
        for m in catalog.GetMatches(mol):
            desc = m.GetDescription() or ""
            fam = (m.GetProp("FilterSet") if m.HasProp("FilterSet") else "") or ""
            key = desc.lower()
            if key in seen:
                continue
            seen.add(key)
            out.append({
                "type": "structural_alert", "severity": "medium",
                "title": f"Structural alert: {desc.replace('_', ' ')}",
                "detail": f"Flagged by the {fam or 'PAINS/Brenk/NIH'} alert set — a substructure "
                          "associated with assay interference, reactivity, or tox liabilities.",
                "action": "Excise or replace the flagged substructure; if it is essential for potency, "
                          "confirm the activity is real with an orthogonal assay (not a PAINS artifact).",
                "tradeoff": "The alerting group may be load-bearing for binding — verify before removing."})
    except Exception as e:
        logger.debug(f"structural alerts skipped: {e}")
    return out[:6]


def analyze(smiles: str, area: str = "", disease: str = "",
            therapeutic_areas: Optional[List[str]] = None) -> Dict:
    """Return a prioritized, structure-aware liability analysis for a molecule."""
    if not _RDKIT:
        return {"available": False, "reason": "RDKit unavailable", "liabilities": []}
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return {"available": False, "reason": "invalid SMILES", "liabilities": []}

    mw   = round(Descriptors.MolWt(mol), 1)
    logp = round(Crippen.MolLogP(mol), 2)
    tpsa = round(rdMolDescriptors.CalcTPSA(mol), 1)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    rtb  = rdMolDescriptors.CalcNumRotatableBonds(mol)
    arom = rdMolDescriptors.CalcNumAromaticRings(mol)
    qed  = round(Descriptors.qed(mol), 3)

    liabilities: List[Dict] = []

    # 1. Reactive / metabolic substructure liabilities (the deep, specific ones)
    for name, smarts, sev, why, action, tradeoff, rxn in _LIABILITY_RULES:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None and mol.HasSubstructMatch(patt):
            liabilities.append({
                "type": "substructure", "severity": sev, "title": name,
                "detail": why, "action": action, "tradeoff": tradeoff,
                # one-click applicable fix (RDKit reaction SMARTS); "" = advice only
                "rxn": rxn, "applicable": bool(rxn)})

    # 2. Structural-alert catalogs (PAINS / Brenk / NIH)
    liabilities.extend(_structural_alerts(mol))

    # 3. Property-level liabilities, framed by mechanism
    liabilities.extend(_limiting_property_advice(mw, logp, tpsa, hbd, hba, rtb, arom))

    # ── CNS-MPO limiting term (only surfaced for CNS/neuro context) ───────────
    cns = None
    is_cns = any(k in (f"{area} {disease} {' '.join(therapeutic_areas or [])}").lower()
                 for k in ("cns", "neuro", "brain", "alzheim", "parkinson", "epileps",
                           "depress", "schizo", "migraine", "psych", "dementia"))
    if is_cns:
        try:
            from services.cns_mpo import cns_mpo
            cns = cns_mpo(smiles=smiles)
            comps = cns.get("components", {}) if cns else {}
            worst = min((c for c in comps.values() if c.get("desirability") is not None),
                        key=lambda c: c["desirability"], default=None)
            worst_name = next((k for k, v in comps.items() if v is worst), "") if worst else ""
            if worst and worst["desirability"] < 0.5 and cns.get("score") is not None:
                liabilities.append({
                    "type": "cns_mpo", "severity": "high" if cns["score"] < 4 else "medium",
                    "property": worst_name,
                    "title": f"CNS-MPO limited by {worst_name} = {worst['value']}",
                    "detail": f"CNS-MPO {cns['score']}/6 ({cns['verdict']}). The single term dragging "
                              f"brain penetration down is {worst_name} (desirability "
                              f"{worst['desirability']:.2f}) — fixing it moves the MPO most.",
                    "action": f"Optimize {worst_name} specifically toward the CNS-desirable range; the "
                              "other MPO terms are already acceptable, so target this one.",
                    "tradeoff": "Balance against the other five MPO terms — over-correcting one can "
                                "push another out of range."})
        except Exception as e:
            logger.debug(f"cns_mpo in advisor skipped: {e}")

    # ── Area-aware developability (which property fails FOR THIS ROUTE) ───────
    develop = None
    try:
        from services import developability
        develop = developability.score(smiles, area=area, therapeutic_areas=therapeutic_areas)
        for p in (develop.get("properties", []) if develop else []):
            if p.get("status") == "out":
                liabilities.append({
                    "type": "developability", "severity": "medium", "property": p.get("key"),
                    "title": f"{p.get('label')} out of range for {develop.get('profile_label','this route')}",
                    "detail": f"{p.get('label')} = {p.get('value')}{p.get('unit','')} is outside the "
                              f"{develop.get('route','') or 'route'} target ({p.get('target')}). "
                              "Developability is route-specific — this is the constraint that matters here.",
                    "action": f"Bring {p.get('label')} into {p.get('target')} for the "
                              f"{develop.get('profile_label','target')} profile.",
                    "tradeoff": "Route-fit changes can trade off against target potency."})
    except Exception as e:
        logger.debug(f"developability in advisor skipped: {e}")

    # De-dupe by title, then rank: severity desc, substructure/alert before property.
    _TYPE_RANK = {"substructure": 3, "structural_alert": 2, "cns_mpo": 2,
                  "developability": 1, "property": 1}
    seen, uniq = set(), []
    for l in liabilities:
        t = l["title"].lower()
        if t in seen:
            continue
        seen.add(t)
        uniq.append(l)
    uniq.sort(key=lambda l: (_SEV_RANK.get(l["severity"], 0),
                             _TYPE_RANK.get(l["type"], 0)), reverse=True)

    limiting = uniq[0] if uniq else None
    if not uniq:
        summary = "No major structural or property liabilities detected — this molecule is a clean " \
                  "starting point; optimize on potency/selectivity rather than ADME-Tox."
    else:
        n_high = sum(1 for l in uniq if l["severity"] == "high")
        lead = limiting["title"]
        summary = (f"{len(uniq)} liability(ies) found ({n_high} high-severity). The one to fix first "
                   f"is: {lead}.")

    return {
        "available": True,
        "properties": {"mw": mw, "logp": logp, "tpsa": tpsa, "hbd": hbd, "hba": hba,
                       "rot_bonds": rtb, "arom_rings": arom, "qed": qed},
        "cns_mpo": cns,
        "developability": {k: develop.get(k) for k in ("profile_label", "route", "score", "level")}
                          if develop and develop.get("available") else None,
        "liabilities": uniq,
        "limiting_liability": limiting,
        "n_high": sum(1 for l in uniq if l["severity"] == "high"),
        "summary": summary,
    }


if __name__ == "__main__":
    import json
    tests = {
        "nitro+ester (tolcapone-like)": "CC(=O)Oc1ccc(C(=O)c2ccc(O)c([N+](=O)[O-])c2)cc1",
        "clean CNS (donepezil)":        "O=C1CC2(CCN(Cc3ccccc3)CC2)Cc2cc(OC)c(OC)cc21",
        "high logP":                    "CCCCCCCCCCCCc1ccc(-c2ccc(-c3ccccc3)cc2)cc1",
    }
    for name, smi in tests.items():
        r = analyze(smi, area="CNS", disease="Alzheimer Disease")
        print(f"\n### {name}")
        print("  ", r["summary"])
        for l in r["liabilities"][:4]:
            print(f"   [{l['severity']}] {l['title']} -> {l['action'][:70]}")
