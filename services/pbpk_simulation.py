"""
Physiologically-Based Pharmacokinetic (PBPK) Simulation Service
================================================================

A genuine perfusion-limited, multi-organ PBPK model (not a cosmetic 2-compartment
fit). Each organ is a real compartment with a mass-balance ODE driven by its
physiological volume and blood flow:

        V_t · dC_t/dt = Q_t · (C_blood − C_t / Kp_t)

with portal first-pass (gut → liver), well-stirred hepatic clearance, and
glomerular renal clearance. ADME inputs are grounded in first principles:

  • fu     — plasma free fraction, DECREASING with lipophilicity
  • Kp     — tissue:plasma partition via Poulin–Theil tissue composition
  • CL_h   — well-stirred hepatic model (bounded by liver blood flow)
  • CL_r   — glomerular filtration of free drug (GFR · fu)
  • F      — oral bioavailability emerges from gut absorption + hepatic first-pass

Efficacy is judged on FREE drug at the target tissue versus the docking-derived
Ki: target occupancy = C_free / (C_free + Ki), plus time above Ki.

Outputs are first-principles ESTIMATES (clearly labelled), not a substitute for a
measured clinical PK study.
"""

import json
import logging
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# ── Measured human PK reference (editable; anchors the model when a drug is present) ──
_PK_REF_FILE = Path(__file__).parent.parent / "data" / "human_pk_reference.json"
_pk_ref_cache: Optional[dict] = None


def _load_pk_reference() -> dict:
    global _pk_ref_cache
    if _pk_ref_cache is None:
        try:
            _pk_ref_cache = json.loads(_PK_REF_FILE.read_text(encoding="utf-8"))
        except Exception as e:
            logger.debug(f"human PK reference unavailable: {e}")
            _pk_ref_cache = {}
    return _pk_ref_cache


def _lookup_human_pk(drug_name: str) -> Optional[dict]:
    """Return measured human PK for a drug (case-insensitive, salt-tolerant), else None."""
    if not drug_name:
        return None
    drugs = _load_pk_reference().get("drugs", {})
    key = drug_name.strip().lower()
    if key in drugs:
        return drugs[key]
    for k, v in drugs.items():       # tolerate salts / suffixes (e.g. "erlotinib hydrochloride")
        if k in key or key in k:
            return v
    return None

try:
    from services.disease_config import get_target_tissues, get_disease_category
    DISEASE_CONFIG_AVAILABLE = True
except Exception:
    DISEASE_CONFIG_AVAILABLE = False
    def get_target_tissues(_): return ["plasma"]
    def get_disease_category(_): return "general"

R_GAS = 0.0019872          # kcal/(mol·K)
BODY_TEMP_K = 310.0        # 37 °C
RT = R_GAS * BODY_TEMP_K   # ≈ 0.616 kcal/mol

# ── Physiology: 70 kg adult (Brown et al. 1997, ICRP). Volumes in L, flows in L/h ──
# Hepatic inflow = hepatic artery (Q_ha) + portal (gut effluent, Q_gut).
_ORGANS = {
    #          V(L)    Q(L/h)
    "liver":  (1.80,  26.0),   # hepatic artery only (portal added via gut routing)
    "gut":    (1.65,  60.0),   # portal flow → liver
    "kidney": (0.31,  74.0),
    "brain":  (1.40,  44.0),
    "muscle": (28.0,  45.0),
    "fat":    (14.0,  16.0),
    "heart":  (0.33,  16.0),
    "lung":   (0.50,   6.0),   # bronchial (tissue) perfusion
    "skin":   (2.60,  18.0),
    "rest":   (14.0,  85.0),
}
_BLOOD_VOLUME = 5.4   # L
_GFR = 7.5            # L/h (≈125 mL/min)

# Poulin–Theil tissue composition: (water, neutral lipid, neutral phospholipid) fractions.
_TISSUE_COMP = {
    "plasma": (0.945, 0.0035, 0.00225),
    "liver":  (0.751, 0.0348, 0.0252),
    "gut":    (0.718, 0.0375, 0.0124),
    "kidney": (0.783, 0.0123, 0.0284),
    "brain":  (0.770, 0.0392, 0.0533),
    "muscle": (0.760, 0.0049, 0.0072),
    "fat":    (0.180, 0.7900, 0.0020),
    "heart":  (0.758, 0.0115, 0.0166),
    "lung":   (0.811, 0.0030, 0.0090),
    "skin":   (0.718, 0.0284, 0.0111),
    "rest":   (0.760, 0.0100, 0.0120),
}

# Disease tissue keyword → PBPK compartment
_TISSUE_ALIAS = {
    "brain": "brain", "cns": "brain", "substantia_nigra": "brain",
    "liver": "liver", "kidney": "kidney", "renal": "kidney",
    "lung": "lung", "respiratory": "lung", "bronchi": "lung", "airway": "lung",
    "gut": "gut", "colon": "gut", "heart": "heart", "vasculature": "heart",
    "skin": "skin", "integument": "skin", "muscle": "muscle",
    "joints": "muscle", "synovium": "muscle", "bone": "rest", "pancreas": "rest",
    "immune": "plasma", "lymphoid": "plasma", "breast": "rest", "tumor": "rest",
}


# Route-of-administration models. k_abs (1/h) and systemic availability f are
# literature-typical ESTIMATES; first-pass metabolism applies only to enteral
# routes. `entry` = the compartment absorbed drug enters. Intranasal carries a
# small direct nose-to-brain fraction that bypasses the blood-brain barrier.
_ROUTE_MODEL: Dict[str, dict] = {
    "iv":          {"label": "Intravenous (IV)",   "entry": "blood", "k_abs": None, "f": 1.00, "first_pass": False},
    "oral":        {"label": "Oral",               "entry": "gut",   "k_abs": "ka", "f": None, "first_pass": True},
    "sublingual":  {"label": "Sublingual",         "entry": "blood", "k_abs": 1.2,  "f": 0.70, "first_pass": False},
    "buccal":      {"label": "Buccal",             "entry": "blood", "k_abs": 0.8,  "f": 0.60, "first_pass": False},
    "sc":          {"label": "Subcutaneous (SC)",  "entry": "blood", "k_abs": 0.4,  "f": 0.90, "first_pass": False},
    "im":          {"label": "Intramuscular (IM)", "entry": "blood", "k_abs": 0.8,  "f": 0.95, "first_pass": False},
    "intranasal":  {"label": "Intranasal",         "entry": "blood", "k_abs": 1.5,  "f": 0.60, "first_pass": False, "nose_to_brain": 0.10},
    "inhalation":  {"label": "Inhalation",         "entry": "blood", "k_abs": 2.0,  "f": 0.60, "first_pass": False},
    "transdermal": {"label": "Transdermal",        "entry": "blood", "k_abs": 0.05, "f": 0.50, "first_pass": False},
    "rectal":      {"label": "Rectal",             "entry": "blood", "k_abs": 0.5,  "f": 0.60, "first_pass": False},
    "ocular":      {"label": "Ophthalmic (local)", "entry": "blood", "k_abs": 0.3,  "f": 0.05, "first_pass": False},
    "topical":     {"label": "Topical (local)",    "entry": "skin",  "k_abs": 0.15, "f": 0.05, "first_pass": False},
}


class PBPKSimulator:
    """Perfusion-limited whole-body PBPK model with free-drug target engagement."""

    def __init__(self, disease_name: str = ""):
        self.disease_name = disease_name
        self.human_params = {"body_weight": 70.0, "blood_volume": _BLOOD_VOLUME}
        self.target_tissue = self._resolve_target_tissue(disease_name)
        logger.info(f"PBPK initialised for {disease_name} (target tissue: {self.target_tissue})")

    # ── Disease / tissue configuration ────────────────────────────────────────
    def _resolve_target_tissue(self, disease_name: str) -> str:
        try:
            tissues = get_target_tissues(disease_name) or []
            for t in tissues:
                comp = _TISSUE_ALIAS.get(t.lower().strip())
                if comp:
                    return comp
        except Exception:
            pass
        return "plasma"

    def set_disease(self, disease_name: str):
        self.disease_name = disease_name
        self.target_tissue = self._resolve_target_tissue(disease_name)

    def get_primary_target_tissue(self) -> str:
        return self.target_tissue

    # ── ADME estimation (first-principles) ────────────────────────────────────
    @staticmethod
    def _free_fraction(logp: float) -> float:
        """Plasma free fraction — decreases with lipophilicity (protein binding rises)."""
        fu = 1.0 / (1.0 + 10.0 ** (0.7 * (logp - 1.5)))
        return float(min(0.99, max(0.01, fu)))

    @staticmethod
    def _kp_poulin_theil(logp: float, fu: float) -> Dict[str, float]:
        """Tissue:plasma partition coefficients via Poulin–Theil composition model."""
        P = 10.0 ** logp
        def num(comp):
            fw, fl, fp = comp
            return P * (fl + 0.3 * fp) + (fw + 0.7 * fp)
        num_p = num(_TISSUE_COMP["plasma"])
        kp = {}
        for t, comp in _TISSUE_COMP.items():
            if t == "plasma":
                continue
            kpu = num(comp) / num_p           # unbound partition
            kp[t] = max(0.03, kpu * fu)        # total tissue:plasma
        return kp

    def _estimate_adme(self, mw: float, logp: float, psa: float, hbd: int,
                       binding_affinity: Optional[float], drug_name: str = "") -> Dict:
        # ── First-principles estimates ────────────────────────────────────────
        fu = self._free_fraction(logp)
        kp = self._kp_poulin_theil(logp, fu)
        q_hep = _ORGANS["liver"][1] + _ORGANS["gut"][1]            # hepatic artery + portal
        clint = 10.0 ** (0.5 * logp + 1.2)
        e_h = (fu * clint) / (q_hep + fu * clint)
        cl_h = e_h * q_hep
        cl_r = _GFR * fu
        ka = min(2.0, max(0.1, 1.5 - psa / 150.0 - hbd * 0.1))
        fa = min(0.99, max(0.10, 1.0 - psa / 250.0 - mw / 4000.0))
        vss_est = _BLOOD_VOLUME + sum(_ORGANS[t][0] * kp[t] for t in _ORGANS)

        anchored = False
        pk_source = "first-principles estimate (no measured PK on record)"

        # ── Anchor to measured human PK when available ────────────────────────
        ref = _lookup_human_pk(drug_name)
        if ref:
            anchored = True
            pk_source = ref.get("source", "measured human PK")
            fu = float(ref.get("fu", fu))
            cl_total = float(ref["cl_l_h"])
            fe = ref.get("fe_renal")
            cl_r = (float(fe) * cl_total) if fe is not None else min(cl_total, _GFR * fu)
            cl_h = max(0.0, cl_total - cl_r)
            # Inverse well-stirred → CLint reproducing the measured hepatic clearance
            clint = (cl_h * q_hep / (fu * (q_hep - cl_h))) if (cl_h < 0.98 * q_hep and fu > 0) else 1e7
            e_h = cl_h / q_hep if q_hep else 0.0
            # Anchor Vd: scale all Kp so Vblood + ΣV·Kp == measured Vss (tissue RATIOS preserved)
            vss_meas = float(ref["vss_l"])
            scale = max(0.02, (vss_meas - _BLOOD_VOLUME) / max(0.1, (vss_est - _BLOOD_VOLUME)))
            kp = {t: v * scale for t, v in kp.items()}
            # Anchor F via fraction absorbed, and ka via measured Tmax
            f_meas = ref.get("f_oral")
            if f_meas is not None:
                fa = min(0.99, max(0.05, float(f_meas) / max(0.05, 1.0 - e_h)))
            if ref.get("tmax_h"):
                ka = min(3.0, max(0.1, 2.5 / float(ref["tmax_h"])))
            vss = vss_meas
            f_oral = float(f_meas) if f_meas is not None else fa * (1.0 - e_h)
        else:
            cl_total = cl_h + cl_r
            vss = vss_est
            f_oral = fa * (1.0 - e_h)

        adme = {
            "fu": round(fu, 4), "clint": round(clint, 2), "e_h": round(e_h, 3),
            "cl_hepatic": round(cl_h, 3), "cl_renal": round(cl_r, 3), "cl": round(cl_total, 3),
            "ka": round(ka, 3), "fa": round(fa, 3), "f": round(f_oral, 3),
            "vss_l": round(vss, 2), "vd": round(vss / self.human_params["body_weight"], 3),
            "ke": round(cl_total / max(vss, 0.1), 4),
            "anchored": anchored, "pk_source": pk_source,
            "measured_t_half": (float(ref["t_half_h"]) if (anchored and ref and ref.get("t_half_h")) else None),
        }
        for t in _ORGANS:
            adme[f"kp_{t}"] = round(kp[t], 3)
        adme["kp"] = kp
        if binding_affinity and binding_affinity < -9.0:
            adme["high_affinity"] = True
        return adme

    # ── ODE integration (perfusion-limited, portal first-pass) ────────────────
    def _integrate(self, dose_mg: float, duration_h: float, route: str, adme: Dict
                   ) -> Tuple[List[float], Dict[str, List[float]]]:
        kp    = adme["kp"]
        cl_r  = adme["cl_renal"]
        clint = adme["clint"]
        fu    = adme["fu"]
        ka    = adme["ka"]
        fa    = adme["fa"]

        tissues = list(_ORGANS.keys())
        m  = len(tissues)
        gi = tissues.index("gut")
        li = tissues.index("liver")
        si = tissues.index("skin")
        V_arr  = np.array([_ORGANS[t][0] for t in tissues])
        Q_arr  = np.array([_ORGANS[t][1] for t in tissues])
        kp_arr = np.array([kp[t] for t in tissues])
        Vb   = _BLOOD_VOLUME
        q_ha = Q_arr[li]
        q_gu = Q_arr[gi]
        q_co = float(Q_arr.sum())

        # Absorption rate + initial conditions per route (literature-typical models)
        rm = _ROUTE_MODEL.get(route, _ROUTE_MODEL["oral"])
        entry = rm["entry"]
        if route == "oral":
            k_abs, f_abs = ka, fa
        elif route == "iv":
            k_abs, f_abs = 0.0, 1.0
        else:
            k_abs, f_abs = float(rm["k_abs"]), float(rm["f"])
        n2b = float(rm.get("nose_to_brain", 0.0))
        bi = tissues.index("brain")

        y0 = np.zeros(m + 2)                       # [Xblood, Xtissue..., depot]
        if route == "iv":
            y0[0] = dose_mg
        else:
            y0[-1] = f_abs * dose_mg               # bioavailable fraction held in depot
            if n2b > 0:                            # intranasal: direct nose-to-brain bolus
                brain_dose = n2b * dose_mg
                y0[1 + bi] = brain_dose
                y0[-1] = max(0.0, y0[-1] - brain_dose)

        def f(t, y):
            Xb = y[0]; depot = y[-1]
            Xt = y[1:1 + m]
            Cb = Xb / Vb
            E  = (Xt / V_arr) / kp_arr             # effluent blood conc per tissue
            dXt = Q_arr * (Cb - E)                 # base perfusion exchange
            absorb = k_abs * depot
            # gut: portal effluent goes to liver; oral absorption enters gut (then first-pass)
            dXt[gi] = Q_arr[gi] * (Cb - E[gi]) + (absorb if route == "oral" else 0.0)
            # liver: hepatic artery + portal in − combined out − intrinsic clearance
            liver_out = (q_ha + q_gu) * E[li]
            dXt[li] = q_ha * Cb + q_gu * E[gi] - liver_out - clint * fu * E[li]
            # local-skin routes (topical): absorption stays in skin
            if entry == "skin":
                dXt[si] = dXt[si] + absorb
            # blood: venous returns (excl. gut/liver) + liver effluent − arterial out − renal CL
            venous_direct = float(np.sum(Q_arr * E) - Q_arr[gi] * E[gi] - Q_arr[li] * E[li])
            dXb = venous_direct + liver_out - q_co * Cb - cl_r * Cb
            # parenteral / nasal / transdermal / inhalation: absorbed drug enters blood, no first-pass
            if entry == "blood" and route != "iv":
                dXb = dXb + absorb
            d_depot = -k_abs * depot
            out = np.empty(m + 2)
            out[0] = dXb; out[1:1 + m] = dXt; out[-1] = d_depot
            return out

        n_eval = max(60, int(duration_h * 20) + 1)
        t_eval = np.linspace(0.0, duration_h, n_eval)
        try:
            from scipy.integrate import solve_ivp
            sol = solve_ivp(f, (0.0, duration_h), y0, method="LSODA",
                            t_eval=t_eval, rtol=1e-6, atol=1e-9, max_step=0.5)
            Y = sol.y
            times = [float(x) for x in sol.t]
        except Exception as e:
            logger.warning(f"solve_ivp failed ({e}); using fixed-step fallback")
            times, Y = self._rk4_fallback(f, y0, duration_h, m)
            times = [float(x) for x in times]

        conc = {"plasma": [float(max(0.0, v / Vb) * 1000.0) for v in Y[0]]}
        for i, t in enumerate(tissues):
            conc[t] = [float(max(0.0, v / V_arr[i]) * 1000.0) for v in Y[1 + i]]
        conc["target"] = conc.get(self.target_tissue, conc["plasma"])
        return times, conc

    @staticmethod
    def _rk4_fallback(f, y0, duration_h, m):
        dt = 0.002
        n = int(duration_h / dt) + 1
        y = y0.copy()
        cols = [[] for _ in range(m + 2)]
        times = []
        for i in range(n):
            t = i * dt
            times.append(t)
            for j in range(m + 2):
                cols[j].append(y[j])
            k1 = f(t, y); k2 = f(t, y + 0.5*dt*k1)
            k3 = f(t, y + 0.5*dt*k2); k4 = f(t, y + dt*k3)
            y = y + dt/6.0 * (k1 + 2*k2 + 2*k3 + k4)
            y = np.maximum(y, 0.0)
        return times, np.array(cols)

    # ── Public API ────────────────────────────────────────────────────────────
    def simulate_drug_exposure(self, drug_name: str, molecular_weight: float = 350.0,
                               logp: float = 2.5, dose_mg: float = 100.0,
                               route: str = "oral", duration_hours: float = 24.0,
                               binding_affinity: Optional[float] = None,
                               params: Optional[Dict] = None) -> Dict:
        p = params or {}
        mw   = float(p.get("mw", molecular_weight))
        logp = float(p.get("logp", logp))
        psa  = float(p.get("psa", 80.0))
        hbd  = int(p.get("hbd", 2))
        route = (route or "oral").lower()
        if route not in _ROUTE_MODEL:
            route = "oral"
        rmodel = _ROUTE_MODEL[route]

        adme = self._estimate_adme(mw, logp, psa, hbd, binding_affinity, drug_name=drug_name)
        # Route-specific systemic bioavailability (oral F emerges from first-pass;
        # parenteral/nasal/etc. use the route model's typical F)
        route_f = 1.0 if route == "iv" else (adme["f"] if route == "oral" else float(rmodel["f"]))
        times, conc = self._integrate(dose_mg, duration_hours, route, adme)

        plasma = conc["plasma"]
        target = conc["target"]
        pk = self._pk_metrics(times, plasma, adme, binding_affinity, mw)

        peak_idx = int(np.argmax(plasma)) if plasma else 0
        plasma_pk = plasma[peak_idx] if plasma else 1.0
        bpr = round((conc["brain"][peak_idx] / plasma_pk), 3) if plasma_pk else 0.0

        # Free-drug target engagement vs docking Ki
        engagement = self._target_engagement(times, target, adme["fu"], binding_affinity, mw)

        step = max(1, int(len(times) / max(1, (duration_hours * 2))))
        organ_ts = {
            "time":   times[::step],
            "plasma": plasma[::step],
            "liver":  conc["liver"][::step],
            "brain":  conc["brain"][::step],
            "kidney": conc["kidney"][::step],
            "muscle": conc["muscle"][::step],
        }

        cl_total = adme["cl"]
        adme_ui = {
            "F_pct":        round(route_f * 100),
            "route_label":  rmodel["label"],
            "t_half":       pk["t_half_hours"],
            "vd_l":         pk["vd_l"],
            "cl_l_h":       pk["clearance_l_h"],
            "hepatic_frac": round(adme["cl_hepatic"] / max(cl_total, 1e-6) * 100),
            "renal_frac":   round(adme["cl_renal"] / max(cl_total, 1e-6) * 100),
            "logp":         round(logp, 2),
            "bbb_ratio":    bpr,
            "fu":           adme["fu"],
            "e_h":          adme["e_h"],
            "anchored":     adme.get("anchored", False),
            "pk_source":    adme.get("pk_source", ""),
        }

        return {
            "drug_name": drug_name, "route": route, "dose_mg": dose_mg,
            "time_hours": times,
            "plasma_concentration_ng_ml": plasma,
            "liver_concentration_ng_ml":  conc["liver"],
            "brain_concentration_ng_ml":  conc["brain"],
            "target_tissue_concentration_ng_ml": target,
            "pk_metrics": pk,
            "cmax_plasma": pk["cmax_ng_ml"], "tmax": pk["tmax_hours"],
            "auc": pk["auc_ng_h_ml"], "half_life": pk["t_half_hours"],
            "brain_plasma_ratio": bpr,
            "organ_timeseries": organ_ts,
            "adme_ui": adme_ui,
            "adme_parameters": adme,
            "target_tissue": self.target_tissue,
            "target_engagement": engagement,
            "model_type": "Perfusion-limited whole-body PBPK (11 organs, portal first-pass)",
            "provenance": (
                (f"CL, Vd, fu and F anchored to measured human PK — {adme['pk_source']}. "
                 "Tissue distribution and free-drug target occupancy from the perfusion-limited model.")
                if adme.get("anchored") else
                ("First-principles PBPK estimate: Poulin–Theil partitioning, well-stirred hepatic "
                 "clearance, glomerular renal clearance. No measured human PK on record for this drug.")),
            "safety_assessment": self._assess_safety(pk, binding_affinity, engagement),
            "success": True,
        }

    def _pk_metrics(self, times, plasma, adme, binding_affinity, mw) -> Dict:
        cmax = max(plasma) if plasma else 0.0
        tmax = times[plasma.index(cmax)] if plasma else 0.0
        auc = 0.0
        for i in range(len(times) - 1):
            auc += (plasma[i] + plasma[i + 1]) / 2.0 * (times[i + 1] - times[i])
        # Measured terminal t½ when anchored; else from the simulated terminal slope
        t_half = (adme.get("measured_t_half")
                  or self._terminal_half_life(times, plasma)
                  or round(0.693 * adme["vss_l"] / max(adme["cl"], 1e-6), 2))
        dt = times[1] - times[0] if len(times) > 1 else 0.0
        thr = 100.0
        time_above = sum(1 for c in plasma if c > thr) * dt
        return {
            "cmax_ng_ml": round(cmax, 2), "tmax_hours": round(tmax, 2),
            "auc_ng_h_ml": round(auc, 2), "t_half_hours": t_half,
            "time_above_threshold_hours": round(time_above, 2),
            "vd_l": adme["vss_l"], "clearance_l_h": round(adme["cl"], 2),
        }

    @staticmethod
    def _terminal_half_life(times, conc) -> Optional[float]:
        """Log-linear regression on the terminal phase → t½."""
        try:
            peak = conc.index(max(conc))
            xs, ys = [], []
            for i in range(peak + 1, len(times)):
                if conc[i] > 1e-6:
                    xs.append(times[i]); ys.append(math.log(conc[i]))
            if len(xs) < 5:
                return None
            tail = max(5, int(len(xs) * 0.5))
            xs, ys = xs[-tail:], ys[-tail:]
            n = len(xs); sx = sum(xs); sy = sum(ys)
            sxx = sum(x*x for x in xs); sxy = sum(x*y for x, y in zip(xs, ys))
            slope = (n*sxy - sx*sy) / (n*sxx - sx*sx)
            if slope >= 0:
                return None
            return round(0.693 / (-slope), 2)
        except Exception:
            return None

    def _target_engagement(self, times, target_conc, fu, binding_affinity, mw) -> Optional[Dict]:
        """Free-drug target occupancy vs docking Ki — the meaningful efficacy readout."""
        if binding_affinity is None or binding_affinity >= 0:
            return None
        ki_M = math.exp(binding_affinity / RT)        # ΔG (kcal/mol, negative) → Ki (M)
        ki_nM = ki_M * 1e9
        occ, free_M = [], []
        for c in target_conc:
            cf_M = fu * c * 1e-6 / mw                  # ng/mL → mol/L, free
            free_M.append(cf_M)
            occ.append(cf_M / (cf_M + ki_M) if (cf_M + ki_M) > 0 else 0.0)
        dt = times[1] - times[0] if len(times) > 1 else 0.0
        time_above_ki = sum(1 for c in free_M if c > ki_M) * dt
        peak_occ = max(occ) if occ else 0.0
        step = max(1, int(len(times) / max(1, (times[-1] * 2 if times else 1))))
        return {
            "ki_nM": round(ki_nM, 3),
            "binding_affinity_kcal_mol": round(binding_affinity, 2),
            "peak_occupancy_pct": round(peak_occ * 100, 1),
            "time_above_ki_hours": round(time_above_ki, 2),
            "target_tissue": self.target_tissue,
            "occupancy_timeseries": {"time": times[::step],
                                     "occupancy_pct": [round(o * 100, 1) for o in occ[::step]]},
            "interpretation": (
                f"Peak free-drug occupancy of {self.target_tissue} target ≈ {round(peak_occ*100)}%; "
                f"free concentration exceeds Ki for {round(time_above_ki,1)} h."),
        }

    def _assess_safety(self, pk, binding_affinity, engagement) -> Dict:
        warnings, margin = [], "Good"
        if pk["cmax_ng_ml"] > 10000:
            margin = "Caution"; warnings.append("High peak plasma concentration")
        if pk["auc_ng_h_ml"] > 200000:
            margin = "Caution"; warnings.append("High total exposure (accumulation risk)")
        if engagement and engagement["peak_occupancy_pct"] < 20:
            warnings.append("Low predicted target occupancy at this dose — efficacy uncertain")
        if not warnings:
            warnings.append("No major exposure flags at the simulated dose")
        return {"safety_margin": margin, "warnings": warnings,
                "therapeutic_window": "Acceptable" if margin == "Good" else "Requires monitoring"}

    # ── Repurposing feasibility (free-drug, occupancy-driven) ──────────────────
    def analyze_repurposing_feasibility(self, drug_name: str, molecular_weight: float,
                                        logp: float, binding_affinity_kcal_mol: float,
                                        target_organ: str, dose_mg: float = 100.0,
                                        route: str = "oral",
                                        original_indication: Optional[str] = None,
                                        new_indication: Optional[str] = None) -> Dict:
        comp = _TISSUE_ALIAS.get((target_organ or "").lower(), self.target_tissue)
        self.target_tissue = comp
        sim = self.simulate_drug_exposure(
            drug_name=drug_name, molecular_weight=molecular_weight, logp=logp,
            dose_mg=dose_mg, route=route, duration_hours=24.0,
            binding_affinity=binding_affinity_kcal_mol)

        eng = sim.get("target_engagement") or {}
        occ = eng.get("peak_occupancy_pct", 0.0)
        t_above = eng.get("time_above_ki_hours", 0.0)
        adme = sim["adme_parameters"]; pk = sim["pk_metrics"]

        # Factor 1: target engagement (50%) — occupancy is the pharmacology
        if occ >= 80:   eng_score, eng_txt = 100, "Saturating target occupancy"
        elif occ >= 50: eng_score, eng_txt = 80,  "Strong target occupancy"
        elif occ >= 20: eng_score, eng_txt = 55,  "Moderate target occupancy"
        elif occ > 0:   eng_score, eng_txt = 30,  "Sub-therapeutic occupancy"
        else:           eng_score, eng_txt = 10,  "Negligible occupancy"

        # Factor 2: duration above Ki (25%)
        cover = min(1.0, t_above / 24.0)
        dur_score = round(30 + 70 * cover)
        # Factor 3: distribution to target tissue (15%)
        kp_t = adme.get(f"kp_{comp}", 1.0)
        dist_score = 100 if kp_t >= 1.5 else 80 if kp_t >= 1.0 else 60 if kp_t >= 0.5 else 35
        # Factor 4: exposure safety (10%)
        safe_score = 100 if pk["cmax_ng_ml"] < 10000 else 60

        feasibility = eng_score*0.5 + dur_score*0.25 + dist_score*0.15 + safe_score*0.10

        if feasibility >= 80:   verdict, color = "Highly Feasible", "success"
        elif feasibility >= 60: verdict, color = "Feasible with Optimization", "info"
        elif feasibility >= 40: verdict, color = "Marginal Feasibility", "warning"
        else:                   verdict, color = "Low Feasibility", "error"

        recs = []
        if occ < 50 and occ > 0:
            recs.append(f"Occupancy is sub-saturating — a higher dose or more frequent schedule "
                        f"would raise free {comp} exposure above Ki for longer.")
        if kp_t < 0.5:
            recs.append(f"Limited partitioning into {comp} (Kp≈{kp_t}); consider a route or "
                        f"formulation that favours local delivery.")
        if not recs:
            recs.append(f"Predicted exposure supports engagement of the {comp} target at {dose_mg} mg.")

        return {
            "success": True, "drug_name": drug_name,
            "original_indication": original_indication or "Unknown",
            "new_indication": new_indication or "Unknown",
            "target_organ": comp, "route": route, "dose_mg": dose_mg,
            "binding_affinity_kcal_mol": round(binding_affinity_kcal_mol, 2),
            "ki_nM": eng.get("ki_nM"),
            "peak_occupancy_pct": occ, "time_above_ki_hours": t_above,
            "feasibility_score": round(feasibility, 1),
            "feasibility_verdict": verdict, "verdict_color": color,
            "feasibility_factors": [
                f"Target occupancy: {eng_txt} ({occ:.0f}%)",
                f"Time above Ki: {t_above:.1f} h of 24 h",
                f"Tissue distribution: Kp({comp}) ≈ {kp_t}",
            ],
            "recommendation": f"{drug_name}: {verdict.lower()} for {new_indication or 'the new indication'} "
                              f"based on predicted free-drug occupancy of the {comp} target.",
            "dose_recommendations": recs,
            "pk_metrics": pk, "adme_parameters": adme,
            "target_engagement": eng,
            "safety_assessment": sim["safety_assessment"],
            "time_hours": sim["time_hours"],
            "plasma_concentration_ng_ml": sim["plasma_concentration_ng_ml"],
            "target_concentration_ng_ml": sim["target_tissue_concentration_ng_ml"],
            "provenance": sim["provenance"],
        }


# Global instance
pbpk_simulator = PBPKSimulator()


# ── Route-of-administration analysis ──────────────────────────────────────────
def _chembl_route_flags(chembl_id: str) -> Dict:
    """Coarse known-administration flags from ChEMBL (oral / parenteral / topical)."""
    out = {"oral": None, "parenteral": None, "topical": None}
    if not chembl_id:
        return out
    try:
        from services import http_client
        r = http_client.get(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json", timeout=8)
        if r and r.ok:
            j = r.json()
            out = {"oral": bool(j.get("oral")), "parenteral": bool(j.get("parenteral")),
                   "topical": bool(j.get("topical"))}
    except Exception as e:
        logger.debug(f"ChEMBL route flags failed for {chembl_id}: {e}")
    return out


def feasible_routes(mw: float, logp: float, psa: float, hbd: int) -> Dict[str, Dict]:
    """Physicochemical feasibility per route (rule-based; clearly an estimate)."""
    drug_like_oral = (mw <= 500 and psa <= 140 and hbd <= 5)
    return {
        "iv":          {"feasible": True, "note": "Parenteral — universal; needs an injectable/soluble form."},
        "im":          {"feasible": True, "note": "Intramuscular depot — broadly feasible."},
        "sc":          {"feasible": True, "note": "Subcutaneous — feasible for small molecules & biologics."},
        "oral":        {"feasible": drug_like_oral,
                        "note": "Oral — drug-like permeability (MW≤500, PSA≤140)." if drug_like_oral
                                else "Oral — permeability-limited (high MW/PSA); bioavailability uncertain."},
        "intranasal":  {"feasible": (mw <= 400 and psa <= 90),
                        "note": "Intranasal — favoured by low MW & moderate polarity; nose-to-brain for CNS."},
        "transdermal": {"feasible": (mw <= 400 and 1 <= logp <= 4 and psa <= 90),
                        "note": "Transdermal — needs MW<400 & logP 1–4 for adequate skin flux."},
        "inhalation":  {"feasible": (mw <= 600),
                        "note": "Inhalation — pulmonary delivery; broad for small molecules."},
        "sublingual":  {"feasible": (mw <= 500 and logp >= 1),
                        "note": "Sublingual — lipophilic, low-dose molecules."},
        "buccal":      {"feasible": (mw <= 500 and logp >= 1),
                        "note": "Buccal — lipophilic, low-dose molecules."},
        "rectal":      {"feasible": True, "note": "Rectal — partial first-pass bypass."},
        "ocular":      {"feasible": True, "note": "Ophthalmic — local eye delivery."},
        "topical":     {"feasible": True, "note": "Topical — local skin delivery."},
    }


def _recommend_route(results: List[Dict], context: Dict) -> Optional[Dict]:
    """Recommend the route for the indication's organ (disease-agnostic)."""
    if not results:
        return None
    pool = [r for r in results if r.get("feasible")] or results
    organ = (context or {}).get("organ", "systemic")

    if organ == "brain":
        # Peak brain delivery is the decisive CNS metric — surfaces intranasal's
        # direct nose-to-brain advantage over routes that only reach it systemically.
        best = max(pool, key=lambda r: r["brain_cmax"])
        why = (f"For a CNS target, highest peak brain exposure (brain Cmax {best['brain_cmax']} ng/mL)"
               + (" — direct nose-to-brain transport bypassing the blood-brain barrier, and non-invasive."
                  if best["route"] == "intranasal" else f" at {best['F_pct']}% bioavailability."))
        return {"route": best["route"], "label": best["label"], "rationale": why}

    # Local-delivery organs: prefer the local route when it ran & is feasible
    local = {"eye": "ocular", "skin": "topical", "lung": "inhalation"}.get(organ)
    if local:
        cand = next((r for r in pool if r["route"] == local), None)
        if cand:
            return {"route": cand["route"], "label": cand["label"],
                    "rationale": f"Local {organ} delivery concentrates drug at the target site with minimal systemic exposure."}

    # Systemic / cardiovascular / GI / hepatic / renal / blood: oral if adequate, else best bioavailability
    oral = next((r for r in pool if r["route"] == "oral" and r.get("feasible")), None)
    if oral and oral["F_pct"] >= 30:
        return {"route": "oral", "label": oral["label"],
                "rationale": f"Oral is feasible and patient-friendly with adequate bioavailability ({oral['F_pct']}%)."}
    best = max(pool, key=lambda r: r["F_pct"])
    return {"route": best["route"], "label": best["label"],
            "rationale": f"Highest systemic bioavailability ({best['F_pct']}%) given limited oral absorption."}


def analyze_routes(drug_name: str, chembl_id: str = "", params: Optional[Dict] = None,
                   dose_mg: float = 100.0, disease_name: str = "",
                   target_organ: str = "", binding_affinity: Optional[float] = None) -> Dict:
    """Run the PBPK across every administration route, with known/feasible flags and a
    route recommended for the indication's organ (disease-agnostic)."""
    params = params or {}
    mw = float(params.get("mw", 350.0)); logp = float(params.get("logp", 2.5))
    psa = float(params.get("psa", 80.0)); hbd = int(params.get("hbd", 2))
    try:
        from services.therapeutic_context import therapeutic_context
        context = therapeutic_context(disease_name)
    except Exception:
        context = {"organ": target_organ or "systemic"}
    if not target_organ:
        target_organ = context.get("organ", "")
    sim = PBPKSimulator(disease_name)
    if target_organ:
        sim.target_tissue = _TISSUE_ALIAS.get(target_organ.lower(), sim.target_tissue)
    known = _chembl_route_flags(chembl_id)
    feas = feasible_routes(mw, logp, psa, hbd)
    order = ["iv", "oral", "sc", "im", "intranasal", "transdermal", "inhalation",
             "sublingual", "ocular", "topical"]
    results = []
    for r in order:
        try:
            res = sim.simulate_drug_exposure(drug_name, dose_mg=dose_mg, route=r,
                                             duration_hours=24.0, binding_affinity=binding_affinity,
                                             params=params)
        except Exception as e:
            logger.debug(f"route {r} sim failed: {e}")
            continue
        bconc = res.get("brain_concentration_ng_ml", []) or []
        times = res.get("time_hours", []) or []
        brain_cmax = max(bconc) if bconc else 0.0
        brain_auc = (sum((bconc[i] + bconc[i + 1]) / 2.0 * (times[i + 1] - times[i])
                         for i in range(len(times) - 1))
                     if (len(bconc) == len(times) and len(bconc) > 1) else 0.0)
        results.append({
            "route": r, "label": _ROUTE_MODEL[r]["label"],
            "F_pct": res["adme_ui"]["F_pct"], "cmax_plasma": round(res["cmax_plasma"], 1),
            "tmax_h": res["tmax"], "half_life_h": res["half_life"],
            "brain_cmax": round(brain_cmax, 1), "brain_auc": round(brain_auc, 1),
            "brain_plasma_ratio": res["brain_plasma_ratio"],
            "first_pass": _ROUTE_MODEL[r]["first_pass"],
            "feasible": feas.get(r, {}).get("feasible"),
            "feasibility_note": feas.get(r, {}).get("note", ""),
        })
    return {
        "drug_name": drug_name, "target_tissue": sim.target_tissue, "dose_mg": dose_mg,
        "known_routes": known, "routes": results, "organ": context.get("organ", "systemic"),
        "recommended": _recommend_route(results, context),
        "disclaimer": ("Per-route absorption rate & bioavailability are literature-typical "
                       "estimates (anchored to measured systemic PK where available); the "
                       "intranasal route includes a heuristic nose-to-brain term. Decision "
                       "support — not a substitute for formal CMC/clinical PK."),
    }
