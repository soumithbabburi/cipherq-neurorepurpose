"""
Within-class vs out-of-class differentiation — the dual comparison required when
filing a drug-repurposing rationale (505(b)(2), patents, regulatory dossiers):

  • Within-class  → "Why THIS drug?"       — benchmark the molecule against the
                     other drugs that hit the SAME target/mechanism, to prove it
                     is the uniquely optimal class member (BBB, clinical maturity).
  • Out-of-class  → "Why THIS mechanism?"   — benchmark the molecule's target
                     against the OTHER druggable targets for the disease, to prove
                     the mechanism is among the best-validated, differentiated ones.

Targets are resolved via the platform's resolve_drug (unions ChEMBL molecule
forms, so salts/parents don't hide the mechanism); peers come from ChEMBL
mechanism on the same target; alternative mechanisms come from Open Targets
disease→target evidence; brain-penetrance from the platform CNS-MPO. Decision
support — not a substitute for formal regulatory/IP analysis.
"""
import logging
from collections import Counter
from typing import Dict, List, Optional, Set

from services import http_client

logger = logging.getLogger(__name__)
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

_mol_cache: Dict[str, Dict] = {}
_target_cache: Dict[str, str] = {}
_compare_cache: Dict[str, Dict] = {}   # full (drug, disease) comparison — repeat loads are instant


def _num(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def _action_word(action: str) -> str:
    a = (action or "").upper()
    if any(k in a for k in ("INHIBIT", "ANTAGONIST", "BLOCKER", "NEGATIVE")):
        return "inhibitor"
    if any(k in a for k in ("AGONIST", "ACTIVATOR", "POSITIVE", "OPENER")):
        return "agonist"
    if "MODULATOR" in a:
        return "modulator"
    return "ligand"


def _mol(chembl_id: str) -> Dict:
    """name / max_phase / smiles for a molecule (cached)."""
    if not chembl_id:
        return {}
    if chembl_id in _mol_cache:
        return _mol_cache[chembl_id]
    out = {"chembl_id": chembl_id, "name": chembl_id, "max_phase": 0, "smiles": "", "first_approval": None}
    try:
        r = http_client.get(f"{CHEMBL}/molecule/{chembl_id}.json", timeout=8)
        if r and r.ok:
            j = r.json()
            out["name"] = j.get("pref_name") or chembl_id
            out["max_phase"] = _num(j.get("max_phase"))
            out["smiles"] = (j.get("molecule_structures") or {}).get("canonical_smiles", "")
            out["first_approval"] = j.get("first_approval")
    except Exception as e:
        logger.debug(f"_mol {chembl_id}: {e}")
    _mol_cache[chembl_id] = out
    return out


def _mpo(smiles: str) -> Dict:
    if not smiles:
        return {}
    try:
        from services.cns_mpo import cns_mpo
        return cns_mpo(smiles=smiles)
    except Exception:
        return {}


def _focal_profile(chembl_id: str, drug_name: str, smiles: str) -> Dict:
    """Robust focal-drug profile (genes union across molecule forms) via resolve_drug."""
    prof = {"chembl_id": chembl_id, "name": drug_name or chembl_id, "genes": [],
            "smiles": smiles, "max_phase": 0, "resolved_id": chembl_id}
    try:
        from services.reverse_repurposing import resolve_drug
        info = resolve_drug(chembl_id) or {}
        prof["genes"] = info.get("targets", []) or []
        prof["smiles"] = smiles or info.get("smiles", "")
        prof["max_phase"] = _num(info.get("max_phase"))
        prof["name"] = drug_name or info.get("name", chembl_id)
        prof["resolved_id"] = info.get("chembl_id", chembl_id)
        prof["known_indications"] = info.get("known_indications", []) or []
    except Exception as e:
        logger.debug(f"_focal_profile {chembl_id}: {e}")
    return prof


def _gene_to_target(gene: str) -> str:
    """Resolve a gene symbol → human single-protein ChEMBL target id (cached)."""
    if not gene:
        return ""
    g = gene.upper()
    if g in _target_cache:
        return _target_cache[g]
    tid = ""
    try:
        r = http_client.get(f"{CHEMBL}/target/search.json",
                            params={"q": gene, "limit": 6, "format": "json"}, timeout=8)
        if r and r.ok:
            for t in r.json().get("targets", []):
                if t.get("organism") == "Homo sapiens" and t.get("target_type") == "SINGLE PROTEIN":
                    tid = t.get("target_chembl_id", ""); break
    except Exception as e:
        logger.debug(f"_gene_to_target {gene}: {e}")
    _target_cache[g] = tid
    return tid


def _peers_on_target(target_chembl_id: str, exclude_ids: Set[str], limit: int = 10):
    """Other molecules with a mechanism on this target → (ids, action_word)."""
    ids, actions = [], []
    try:
        r = http_client.get(f"{CHEMBL}/mechanism.json",
            params={"target_chembl_id": target_chembl_id, "limit": 80, "format": "json"}, timeout=10)
        if r and r.ok:
            for m in r.json().get("mechanisms", []):
                if m.get("action_type"):
                    actions.append(m["action_type"])
                mid = m.get("molecule_chembl_id")
                if mid and mid not in exclude_ids and mid not in ids:
                    ids.append(mid)
    except Exception as e:
        logger.debug(f"_peers_on_target {target_chembl_id}: {e}")
    word = _action_word(actions[0]) if not actions else _action_word(Counter(actions).most_common(1)[0][0])
    return ids[:limit], word


def _representative_drug(gene: str) -> Optional[Dict]:
    """Best-effort highest-phase known drug for a target gene (an alternative mechanism)."""
    tid = _gene_to_target(gene)
    if not tid:
        return None
    try:
        rm = http_client.get(f"{CHEMBL}/mechanism.json",
                            params={"target_chembl_id": tid, "limit": 25, "format": "json"}, timeout=8)
        best = None
        if rm and rm.ok:
            for m in rm.json().get("mechanisms", []):
                mid = m.get("molecule_chembl_id")
                if not mid:
                    continue
                mol = _mol(mid)
                if best is None or (mol.get("max_phase") or 0) > (best.get("max_phase") or 0):
                    best = mol
        if best:
            return {"name": best.get("name"), "chembl_id": best.get("chembl_id"), "max_phase": best.get("max_phase", 0)}
    except Exception as e:
        logger.debug(f"_representative_drug {gene}: {e}")
    return None


OPENFDA_EVENT = "https://api.fda.gov/drug/event.json"


def _half_life(name: str, smiles: str) -> Optional[float]:
    """Predicted terminal t½ (h) — measured-anchored where available, else a PBPK
    estimate from physicochemistry. The key duration-of-action axis."""
    if not smiles:
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        m = Chem.MolFromSmiles(smiles)
        if not m:
            return None
        mw, logp = Descriptors.MolWt(m), Descriptors.MolLogP(m)
        psa, hbd = rdMolDescriptors.CalcTPSA(m), rdMolDescriptors.CalcNumHBD(m)
        from services.pbpk_simulation import PBPKSimulator
        adme = PBPKSimulator()._estimate_adme(mw, logp, psa, hbd, None, drug_name=name)
        tm = adme.get("measured_t_half")
        return round(float(tm), 1) if tm else round(0.693 * adme["vss_l"] / max(adme["cl"], 1e-6), 1)
    except Exception as e:
        logger.debug(f"_half_life {name}: {e}")
        return None


def _n_targets(chembl_id: str) -> Optional[int]:
    """Selectivity proxy: # of distinct mechanism targets (fewer = more selective)."""
    try:
        r = http_client.get(f"{CHEMBL}/mechanism.json",
                            params={"molecule_chembl_id": chembl_id, "limit": 50, "format": "json"}, timeout=6)
        if r and r.ok:
            tids = {m.get("target_chembl_id") for m in r.json().get("mechanisms", []) if m.get("target_chembl_id")}
            return len(tids) or None
    except Exception as e:
        logger.debug(f"_n_targets {chembl_id}: {e}")
    return None


def _serious_pct(name: str) -> Optional[int]:
    """Safety signal: % of FAERS reports flagged serious (higher = worse)."""
    if not name:
        return None
    try:
        r = http_client.get(OPENFDA_EVENT,
            params={"search": f'patient.drug.openfda.generic_name:"{name.lower()}"', "count": "serious"}, timeout=8)
        if r and r.ok:
            res = {str(x.get("term")): x.get("count", 0) for x in r.json().get("results", [])}
            s, ns = res.get("1", 0), res.get("2", 0)
            if s + ns > 0:
                return round(100 * s / (s + ns))
    except Exception as e:
        logger.debug(f"_serious_pct {name}: {e}")
    return None


def _offtarget_ae(target_name: str) -> str:
    """Known selectivity-driven adverse effects for common off-targets."""
    n = (target_name or "").lower()
    if "phosphodiesterase 6" in n or "pde6" in n:
        return "PDE6 cross-reactivity is linked to transient visual disturbances (blue-tinged vision)."
    if "phosphodiesterase 11" in n or "pde11" in n:
        return "PDE11 cross-reactivity is associated with myalgia / back pain."
    if "herg" in n or "potassium channel" in n:
        return "hERG activity is a cardiac QT-prolongation liability."
    if "phosphodiesterase 1" in n:
        return "PDE1 cross-reactivity may affect cardiac/vascular tissue."
    return ""


def _selectivity_profile(chembl_id: str, primary_gene: str) -> Dict:
    """Quantitative selectivity: best measured potency on the primary target vs the
    strongest off-target (ChEMBL pChEMBL). Fold = 10^(Δ pChEMBL)."""
    out = {"available": False}
    primary_tid = _gene_to_target(primary_gene)
    try:
        r = http_client.get(f"{CHEMBL}/activity.json",
            params={"molecule_chembl_id": chembl_id, "pchembl_value__isnull": "false",
                    "limit": 1000, "format": "json"}, timeout=15)
        if not (r and r.ok):
            return out
        # Group by target NAME (the same protein can have several ChEMBL target ids)
        by_name: Dict[str, float] = {}
        primary_name = None
        for a in r.json().get("activities", []):
            tid, pv = a.get("target_chembl_id"), a.get("pchembl_value")
            nm = (a.get("target_pref_name") or "").strip()
            if not nm or pv is None:
                continue
            try:
                pv = float(pv)
            except (TypeError, ValueError):
                continue
            by_name[nm] = max(by_name.get(nm, 0.0), pv)
            if tid and tid == primary_tid:
                primary_name = nm
        if not by_name:
            return out
        if not primary_name:
            primary_name = max(by_name, key=by_name.get)         # most-potent target = primary
        prim_p = by_name[primary_name]
        offs = [(nm, p) for nm, p in by_name.items() if nm != primary_name]
        out = {"available": True, "primary": primary_name,
               "primary_pchembl": round(prim_p, 1), "n_targets_measured": len(by_name)}
        if offs:
            off_name, off_p = max(offs, key=lambda kv: kv[1])
            fold = 10 ** (prim_p - off_p) if off_p > 0 else None
            out.update({"top_offtarget": off_name, "offtarget_pchembl": round(off_p, 1),
                        "selectivity_fold": round(fold) if (fold and fold >= 1) else None,
                        "ae_note": _offtarget_ae(off_name)})
        else:
            out["selectivity_fold"] = None
            out["note"] = "Only the primary target has measured potency on record (clean selectivity profile)."
    except Exception as e:
        logger.debug(f"_selectivity_profile {chembl_id}: {e}")
    return out


def _bbb_evidence(name: str, smiles: str, cns_relevance: str) -> Dict:
    """Weight-of-evidence BBB justification: CNS-MPO + predicted brain:plasma Kp +
    a P-gp efflux heuristic. No unsupported 'it crosses' claim."""
    mpo = _mpo(smiles)
    score = mpo.get("score")
    if score is None:
        return {"available": False}
    kp_brain, efflux = None, "unknown"
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        m = Chem.MolFromSmiles(smiles)
        if m:
            mw, logp = Descriptors.MolWt(m), Descriptors.MolLogP(m)
            psa, hbd = rdMolDescriptors.CalcTPSA(m), rdMolDescriptors.CalcNumHBD(m)
            hba = rdMolDescriptors.CalcNumHBA(m)
            from services.pbpk_simulation import PBPKSimulator
            adme = PBPKSimulator()._estimate_adme(mw, logp, psa, hbd, None, drug_name=name)
            kp_brain = round((adme.get("kp") or {}).get("brain", 0.0), 2)
            # P-gp substrates are typically larger & more polar (Hitchcock/efflux rules of thumb)
            efflux = "high" if (psa >= 90 or (mw >= 400 and hba >= 8)) else ("moderate" if psa >= 70 else "low")
    except Exception as e:
        logger.debug(f"_bbb_evidence {name}: {e}")
    if cns_relevance == "required":
        if score >= 4 and efflux in ("low", "moderate"):
            stmt = (f"Physicochemically BBB-capable (CNS-MPO {score}/6, predicted brain:plasma Kp≈{kp_brain}) "
                    f"with {efflux} predicted P-gp efflux — brain exposure is plausible but should be confirmed "
                    f"with measured CNS PK for this indication.")
        elif score >= 4:
            stmt = (f"Physicochemically BBB-capable (CNS-MPO {score}/6) but {efflux} predicted P-gp efflux "
                    f"limits brain exposure (Kp≈{kp_brain}) — central exposure must be justified with measured data.")
        else:
            stmt = (f"Limited BBB permeability (CNS-MPO {score}/6, Kp≈{kp_brain}) — unlikely to reach the CNS target "
                    f"without a delivery strategy (e.g., intranasal).")
    else:
        lim = "limits" if score < 4 else "implies only modest"
        stmt = (f"Predicted brain:plasma Kp≈{kp_brain} (CNS-MPO {score}/6, {efflux} efflux). For this peripheral "
                f"indication that {lim} central exposure — relevant to the central side-effect profile, not efficacy.")
    return {"available": True, "cns_mpo": score, "kp_brain": kp_brain,
            "efflux_risk": efflux, "statement": stmt}


def _pk_fit(half_life: Optional[float], context: Dict, peers: List[Dict]) -> Dict:
    """Frame the focal half-life against the indication — with the CORRECT dosing
    implication (short t½ → MORE frequent dosing, not less)."""
    if half_life is None:
        return {"statement": ""}
    organ = context.get("organ", "systemic")
    disease = (context.get("disease") or "").lower()
    longer = [p for p in peers if p.get("half_life_h") and p["half_life_h"] > 1.3 * half_life]
    longer_name = longer[0]["name"] if longer else None

    # Mathematically correct dosing implication of the half-life
    if half_life < 4:
        dosing = ("entails multiple-daily (or controlled-release) dosing; its rapid clearance is an asset for "
                  "tight titration and low systemic accumulation, especially in patients with reduced organ clearance")
    elif half_life < 12:
        dosing = "is consistent with roughly twice-daily dosing"
    else:
        dosing = "permits convenient once-daily (or less-frequent) dosing, at the cost of slower washout / accumulation"

    if "erectile" in disease:
        stmt = (f"Predicted t½ {half_life} h {dosing}. A short t½ suits on-demand ED use"
                + (f"; longer-acting members (e.g., {longer_name}) instead enable daily dosing." if longer_name else "."))
    elif organ == "brain":
        stmt = (f"Predicted t½ {half_life} h {dosing}. For a chronic CNS indication, sustained CNS exposure and BBB "
                f"penetration matter more than peak — see BBB evidence.")
    else:
        stmt = f"Predicted t½ {half_life} h {dosing}."
        if longer_name:
            stmt += f" Longer-acting peers (e.g., {longer_name}) offer less-frequent dosing for a chronic indication."
    return {"half_life_h": half_life, "statement": stmt}


def _within_class(prof: Dict, context: Dict) -> Dict:
    """Compare the focal drug against other drugs hitting its primary target,
    differentiating on the axis that matters for THIS indication's organ."""
    genes = prof.get("genes", [])
    if not genes:
        return {"available": False,
                "reason": "No protein target on record — within-class peers cannot be defined by mechanism."}
    gene = genes[0]
    target_id = _gene_to_target(gene)
    if not target_id:
        return {"available": False, "reason": f"Could not resolve a ChEMBL target for {gene}."}

    exclude = {prof["chembl_id"], prof.get("resolved_id", prof["chembl_id"])}
    peer_ids, word = _peers_on_target(target_id, exclude, limit=14)

    fmpo = _mpo(prof.get("smiles", ""))
    focal = {"chembl_id": prof["chembl_id"], "name": prof["name"], "max_phase": prof.get("max_phase", 0),
             "cns_mpo": fmpo.get("score"), "cns_druggable": fmpo.get("cns_druggable"),
             "smiles": prof.get("smiles", ""), "focal": True}
    members = [focal]
    focal_lower = (prof["name"] or "").lower()
    seen = {focal_lower}
    for pid in peer_ids:
        m = _mol(pid)
        nm = (m.get("name") or "").lower()
        # dedupe by name AND drop salt/parent forms of the focal drug itself
        if not nm or nm in seen or (focal_lower and (focal_lower in nm or nm in focal_lower)):
            continue
        seen.add(nm)
        pm = _mpo(m.get("smiles", ""))
        members.append({"chembl_id": pid, "name": m.get("name", pid), "max_phase": m.get("max_phase", 0),
                        "cns_mpo": pm.get("score"), "cns_druggable": pm.get("cns_druggable"),
                        "smiles": m.get("smiles", ""), "focal": False})
        if len(members) >= 9:
            break

    peers = sorted([x for x in members if not x["focal"]], key=lambda x: -(x.get("max_phase") or 0))
    members = [focal] + peers

    # ── Computed differentiation axes (half-life always; selectivity & safety
    #    for peripheral indications, where they actually decide the winner) ─────
    peripheral = (context.get("cns_relevance") != "required")
    for mem in members:
        mem["half_life_h"] = _half_life(mem["name"], mem.get("smiles", ""))
        if peripheral:
            mem["n_targets"] = _n_targets(mem["chembl_id"])
            mem["serious_pct"] = _serious_pct(mem["name"])
            # focal mechanism may be under a salt id → fall back to the resolved gene count
            if mem["focal"] and mem.get("n_targets") is None:
                mem["n_targets"] = len(prof.get("genes", [])) or None
        mem.pop("smiles", None)   # don't ship SMILES to the client

    # ── Auto "why this drug" justification — axis depends on the indication ────
    n = len(members)
    bp = [x for x in members if x.get("cns_druggable") is True]
    approved = [x for x in members if (x.get("max_phase") or 0) >= 4]
    rel = context.get("cns_relevance", "neutral")   # required | neutral | liability
    fmpo = focal.get("cns_mpo")
    bits: List[str] = []
    if rel == "required":                            # CNS indication — BBB is the differentiator
        if focal.get("cns_druggable") is True:
            if len(bp) == 1:
                bits.append(f"the ONLY brain-penetrant molecule (CNS-MPO {fmpo}/6) among {n} {gene} {word}s — decisive for a CNS indication")
            else:
                bits.append(f"one of {len(bp)} of {n} brain-penetrant {gene} {word}s (CNS-MPO {fmpo}/6)")
        elif fmpo is not None:
            bits.append(f"physicochemically borderline for the brain (CNS-MPO {fmpo}/6) versus its {gene} {word} peers — a liability for a CNS target")
    elif rel == "liability" and fmpo is not None:    # peripheral indication — CNS exposure is a risk, not a benefit
        if focal.get("cns_druggable") is False:
            bits.append(f"peripherally restricted (CNS-MPO {fmpo}/6 → low CNS exposure) — minimal central side-effect risk for this peripheral indication")
        else:
            bits.append(f"CNS-exposed (CNS-MPO {fmpo}/6) — central side effects to monitor for this peripheral indication")
    fp = focal.get("max_phase") or 0
    if fp >= 4:
        bits.append("an approved drug — a de-risked safety package available for 505(b)(2) citation")
    elif approved:
        bits.append(f"clinically at phase {fp}; {len(approved)} class member(s) are approved, providing reference safety data")
    if rel != "required":
        fhl = focal.get("half_life_h")
        peer_hls = sorted(m["half_life_h"] for m in members if not m["focal"] and m.get("half_life_h"))
        if fhl and peer_hls:
            med = peer_hls[len(peer_hls) // 2]
            if fhl >= 1.5 * med:
                bits.insert(0, f"longer-acting (predicted t½ {fhl} h vs class median {med:g} h) — a duration-of-action advantage")
            elif fhl <= 0.67 * med:
                bits.append(f"shorter-acting (predicted t½ {fhl} h vs class median {med:g} h)")
            else:
                bits.append(f"comparable duration (predicted t½ {fhl} h, class median {med:g} h)")
        else:
            bits.append("the decisive within-class axes here are duration of action, target selectivity and safety")
    uniqueness = (f"{focal['name']} is " + "; ".join(bits) + "."
                  if bits else
                  f"{focal['name']} shares the {gene} {word} class — differentiation needs additional axes (potency, selectivity, IP).")

    # ── Structured "why this drug?" justification (4 filing-grade axes) ───────
    from datetime import date
    fa = _mol(prof["chembl_id"]).get("first_approval")
    try:
        tenure = date.today().year - int(fa) if fa else None
    except (TypeError, ValueError):
        tenure = None
    focal_serious = focal.get("serious_pct")
    if focal_serious is None:
        focal_serious = _serious_pct(prof["name"])
    peer_serious = sorted(m["serious_pct"] for m in members if not m["focal"] and m.get("serious_pct") is not None)
    class_med_serious = peer_serious[len(peer_serious) // 2] if peer_serious else None

    sbits = []
    if tenure:
        sbits.append(f"~{tenure} years of post-market safety")
    if focal_serious is not None and class_med_serious is not None:
        cmp = "below" if focal_serious < class_med_serious else "above" if focal_serious > class_med_serious else "at"
        sbits.append(f"a serious-AE fraction of {focal_serious}% ({cmp} the class median {class_med_serious}%)")
    elif focal_serious is not None:
        sbits.append(f"a serious-AE fraction of {focal_serious}% (FAERS)")
    better = (tenure and class_med_serious is not None and focal_serious is not None and focal_serious <= class_med_serious)
    safety = {
        "first_approval": fa, "tenure_years": tenure, "serious_pct": focal_serious,
        "class_median_serious": class_med_serious,
        "statement": ((f"{focal['name']} carries " + "; ".join(sbits) + " — "
                       + ("the most de-risked safety package in the class for 505(b)(2) citation."
                          if better else "a substantial, citable post-market safety record.")) if sbits else ""),
    }

    justification = {
        "selectivity": _selectivity_profile(prof["chembl_id"], gene),
        "bbb": _bbb_evidence(prof["name"], prof.get("smiles", ""), rel),
        "safety": safety,
        "pk": _pk_fit(focal.get("half_life_h"), context, [m for m in members if not m["focal"]]),
    }
    parts = [uniqueness]
    sel = justification["selectivity"]
    if sel.get("available") and sel.get("selectivity_fold"):
        s = f"It is ~{sel['selectivity_fold']}× selective for {sel['primary']} over {sel.get('top_offtarget', 'off-targets')}"
        if sel.get("ae_note"):
            s += f" ({sel['ae_note']})"
        parts.append(s + ".")
    elif sel.get("note"):
        parts.append(sel["note"])
    for ax in ("bbb", "pk", "safety"):
        st = justification[ax].get("statement")
        if st:
            parts.append(st)
    justification["summary"] = " ".join(parts)

    return {"available": True,
            "class_definition": f"{gene} {word}s (shared primary molecular target, ChEMBL)",
            "target_gene": gene, "n_members": n, "members": members, "uniqueness": uniqueness,
            "justification": justification,
            "cns_relevance": rel,
            "cns_mpo_label": context.get("cns_mpo_label", "CNS-MPO"),
            "cns_mpo_higher_is_better": context.get("cns_mpo_higher_is_better", True)}


def _out_of_class(prof: Dict, disease: str, efo_id: str = "") -> Dict:
    """Compare the focal mechanism against the disease's other druggable targets."""
    focal_genes = {g.upper() for g in prof.get("genes", [])}
    try:
        from services.biocypher_graph import _ot_disease_targets
    except Exception as e:
        logger.debug(f"out-of-class import: {e}")
        return {"available": False, "reason": "Open Targets module unavailable."}

    if not efo_id:
        try:
            from services.disease_ontology import resolve_disease
            efo_id = (resolve_disease(disease) or {}).get("disease_id", "")
        except Exception:
            efo_id = ""
    alts = []
    if efo_id:
        try:
            alts = _ot_disease_targets(efo_id, n=25)
        except Exception as e:
            logger.debug(f"ot targets: {e}")
    if not alts:
        return {"available": False,
                "reason": "No Open Targets disease-target associations found — alternative mechanisms cannot be ranked."}

    ranked = sorted(alts, key=lambda a: -a.get("score", 0))
    focal_rank, focal_score, focal_gene = None, None, ""
    for i, a in enumerate(ranked):
        if a.get("gene_symbol", "").upper() in focal_genes:
            focal_rank, focal_score, focal_gene = i + 1, a.get("score", 0), a.get("gene_symbol", "")
            break

    alternatives = []
    for a in ranked[:6]:
        g = a.get("gene_symbol", "")
        is_focal = g.upper() in focal_genes
        alternatives.append({
            "gene": g, "score": round(a.get("score", 0), 3), "is_focal_mechanism": is_focal,
            "representative_drug": (None if is_focal else _representative_drug(g)),
        })

    other_top = [a["gene"] for a in alternatives if not a["is_focal_mechanism"]][:3]
    if focal_rank:
        tier = "top-tier" if focal_rank <= 3 else "well-supported" if focal_rank <= 10 else "secondary"
        validity = (f"{prof['name']}'s target {focal_gene} ranks #{focal_rank} of {len(ranked)} druggable targets "
                    f"associated with {disease} (Open Targets score {focal_score:.2f}) — a {tier} mechanism"
                    + (f", differentiated from alternatives such as {', '.join(other_top)}." if other_top else "."))
    else:
        validity = (f"{prof['name']}'s target(s) {', '.join(sorted(focal_genes)) or '(none annotated)'} are NOT among the "
                    f"top {len(ranked)} Open Targets-associated targets for {disease} — consistent with a symptomatic / "
                    f"functional mechanism rather than a genetic-risk target. The rationale should lean on clinical and "
                    f"pathway evidence, a point reviewers will probe.")

    return {"available": True, "disease": disease, "efo_id": efo_id,
            "focal_genes": sorted(focal_genes), "focal_rank": focal_rank,
            "focal_score": round(focal_score, 3) if focal_score is not None else None,
            "total_targets": len(ranked), "alternatives": alternatives, "validity": validity}


def _novelty(prof: Dict, disease: str) -> Dict:
    """Phase-aware novelty via the single source of truth (regulatory_verdict).
    A claim only FAILS novelty when the molecule is APPROVED (phase 4) for this
    indication — not when it was merely studied at a low phase or co-mentioned in
    a drug-interaction study. Returns the 3-state verdict (approved_here /
    investigated_here / prior_signal / novel)."""
    try:
        from services.regulatory_verdict import novelty as _nv
        nv = _nv(prof.get("known_indications", []), disease)
        return {
            "is_known_indication": nv["is_known_indication"],
            "state":   nv["state"],
            "matched": nv["matched"],
            "best_phase": nv["best_phase"],
            "note":    nv["note"],
            "patentability": nv.get("patentability", ""),
        }
    except Exception as e:
        logger.debug(f"_novelty delegation failed: {e}")
        return {"is_known_indication": False, "state": "unknown", "matched": None,
                "note": f"Novelty could not be assessed for {prof.get('name','this drug')} in \"{disease}\"."}


def compare_classes(chembl_id: str, disease: str, drug_name: str = "",
                    smiles: str = "", efo_id: str = "") -> Dict:
    """Full dual comparison for a (drug, disease) repurposing hypothesis."""
    key = f"{chembl_id}|{(disease or '').strip().lower()}"
    if key in _compare_cache:
        return _compare_cache[key]
    prof = _focal_profile(chembl_id, drug_name, smiles)
    try:
        from services.therapeutic_context import therapeutic_context
        context = therapeutic_context(disease)
    except Exception:
        context = {"cns_relevance": "neutral", "organ": "systemic",
                   "cns_mpo_label": "CNS-MPO", "cns_mpo_higher_is_better": True}
    result = {
        "chembl_id": chembl_id, "drug_name": prof["name"], "disease": disease,
        "focal_targets": prof.get("genes", []), "context": context,
        "novelty": _novelty(prof, disease),
        "within_class": _within_class(prof, context),
        "out_of_class": _out_of_class(prof, disease, efo_id),
        "disclaimer": ("Within-class = drugs sharing the primary ChEMBL mechanism target; out-of-class = "
                       "Open Targets disease-associated targets. Decision-support for a differentiation "
                       "argument — not a substitute for formal regulatory or freedom-to-operate analysis."),
    }
    _compare_cache[key] = result
    return result
