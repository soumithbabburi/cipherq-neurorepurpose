"""
PBPK (Physiologically-Based Pharmacokinetic) Simulation Service
Human simulation after molecular docking for drug concentration-time predictions
Supports disease-specific pharmacokinetic parameters
"""

import logging
import math
from typing import Dict, List, Tuple, Optional
import json

logger = logging.getLogger(__name__)

# Import disease configuration
try:
    from services.disease_config import (
        get_disease_profile, get_pbpk_adjustments, get_target_tissues,
        get_disease_category
    )
    DISEASE_CONFIG_AVAILABLE = True
except ImportError:
    DISEASE_CONFIG_AVAILABLE = False
    logger.warning("Disease configuration not available for PBPK - using defaults")

class PBPKSimulator:
    """
    Compartmental PBPK model for human drug simulation
    Predicts drug concentration-time curves after docking
    Supports disease-specific tissue targeting and pharmacokinetics
    """
    
    def __init__(self, disease_name: str = "Alzheimer's Disease"):
        self.disease_name = disease_name
        self.disease_category = 'neurological'
        self.target_tissues = ['brain']
        self.pbpk_adjustments = {}
        
        # Default human physiological parameters (70 kg adult)
        self.human_params = {
            'body_weight': 70.0,  # kg
            'blood_volume': 5.0,  # L
            'cardiac_output': 6.0,  # L/min
            'liver_volume': 1.8,  # L
            'kidney_volume': 0.31,  # L
            'gut_volume': 1.65,  # L
            'fat_volume': 14.0,  # L
            'muscle_volume': 28.0,  # L
            'brain_volume': 1.4,  # L
            'lung_volume': 1.2,  # L
            'heart_volume': 0.3,  # L
            'pancreas_volume': 0.1,  # L
            'liver_blood_flow': 1.45,  # L/min (24% of cardiac output)
            'kidney_blood_flow': 1.24,  # L/min (21% of cardiac output)
            'gut_blood_flow': 1.1,  # L/min
            'muscle_blood_flow': 0.75,  # L/min
            'fat_blood_flow': 0.26,  # L/min
            'brain_blood_flow': 0.74,  # L/min
            'lung_blood_flow': 6.0,  # L/min (full cardiac output)
            'heart_blood_flow': 0.24,  # L/min
            'pancreas_blood_flow': 0.06  # L/min
        }
        
        # Configure for disease
        self._configure_for_disease(disease_name)
        
        logger.info(f"PBPK Simulator initialized for {disease_name}")
    
    def _configure_for_disease(self, disease_name: str):
        """Configure PBPK simulator for specific disease"""
        self.disease_name = disease_name
        
        if DISEASE_CONFIG_AVAILABLE:
            profile = get_disease_profile(disease_name)
            self.disease_category = profile.category
            self.target_tissues = profile.target_tissues
            self.pbpk_adjustments = profile.pbpk_adjustments
            logger.info(f"✅ PBPK configured for {disease_name} - targets: {self.target_tissues}")
        else:
            # Default settings
            self.disease_category = 'neurological'
            self.target_tissues = ['brain']
            self.pbpk_adjustments = {'kp_brain': 1.5}
    
    def set_disease(self, disease_name: str):
        """Update disease configuration"""
        self._configure_for_disease(disease_name)
    
    def get_primary_target_tissue(self) -> str:
        """Get the primary target tissue for the current disease"""
        if self.target_tissues:
            return self.target_tissues[0]
        return 'plasma'
    
    def simulate_drug_exposure(
        self,
        drug_name: str,
        molecular_weight: float,
        logp: float,
        dose_mg: float = 100.0,
        route: str = 'oral',
        duration_hours: float = 24.0,
        binding_affinity: Optional[float] = None
    ) -> Dict:
        """
        Simulate drug concentration-time profile in human
        
        Args:
            drug_name: Name of the drug
            molecular_weight: Molecular weight (g/mol)
            logp: Lipophilicity (LogP value)
            dose_mg: Dose in milligrams
            route: Administration route (oral, IV)
            duration_hours: Simulation duration
            binding_affinity: Binding affinity from docking (kcal/mol)
        
        Returns:
            Dictionary with simulation results
        """
        logger.info(f"Starting PBPK simulation for {drug_name}")
        
        # Estimate ADME parameters from molecular properties
        adme_params = self._estimate_adme_parameters(molecular_weight, logp, binding_affinity)
        
        # Set up compartmental model
        if route == 'oral':
            time_points, concentrations = self._simulate_oral_dosing(
                dose_mg, duration_hours, adme_params
            )
        else:  # IV
            time_points, concentrations = self._simulate_iv_dosing(
                dose_mg, duration_hours, adme_params
            )
        
        # Calculate PK parameters
        pk_metrics = self._calculate_pk_metrics(time_points, concentrations, adme_params)
        
        # Convert to lists if they are numpy arrays
        import numpy as np
        time_list = time_points.tolist() if isinstance(time_points, np.ndarray) else time_points
        plasma_list = concentrations['plasma'].tolist() if isinstance(concentrations['plasma'], np.ndarray) else concentrations['plasma']
        liver_list = concentrations['liver'].tolist() if isinstance(concentrations['liver'], np.ndarray) else concentrations['liver']
        brain_list = concentrations['brain'].tolist() if isinstance(concentrations['brain'], np.ndarray) else concentrations['brain']
        target_list = concentrations['target'].tolist() if isinstance(concentrations['target'], np.ndarray) else concentrations['target']
        
        return {
            'drug_name': drug_name,
            'route': route,
            'dose_mg': dose_mg,
            'time_hours': time_list,
            'plasma_concentration_ng_ml': plasma_list,
            'liver_concentration_ng_ml': liver_list,
            'brain_concentration_ng_ml': brain_list,
            'target_tissue_concentration_ng_ml': target_list,
            'pk_metrics': pk_metrics,
            'adme_parameters': adme_params,
            'safety_assessment': self._assess_safety(pk_metrics, binding_affinity),
            'success': True
        }
    
    def _estimate_adme_parameters(
        self,
        molecular_weight: float,
        logp: float,
        binding_affinity: Optional[float]
    ) -> Dict:
        """Estimate ADME parameters from molecular properties with disease-specific adjustments"""
        
        # Absorption rate constant (h^-1)
        ka = 0.5 + (logp / 10.0)  # Lipophilic drugs absorb faster
        ka = max(0.1, min(2.0, ka))
        
        # Apply disease-specific gut absorption adjustment
        if 'gut_absorption' in self.pbpk_adjustments:
            ka *= self.pbpk_adjustments['gut_absorption']
        
        # Volume of distribution (L/kg)
        vd = 0.7 + (logp * 0.3)  # Lipophilic drugs distribute more
        vd = max(0.2, min(4.0, vd))
        
        # Clearance (L/h/kg)
        cl = 0.05 + (1.0 / (molecular_weight / 100.0))  # Smaller molecules clear faster
        cl = max(0.01, min(0.5, cl))
        
        # Apply renal clearance adjustment for kidney diseases
        if 'renal_clearance' in self.pbpk_adjustments:
            cl *= self.pbpk_adjustments['renal_clearance']
        
        # Bioavailability (fraction)
        f = 0.7 - (molecular_weight / 2000.0)  # Large molecules have lower oral bioavailability
        f = max(0.1, min(1.0, f))
        
        # Protein binding (fraction)
        fu = 0.3 + (logp / 15.0)  # Lipophilic drugs bind more to proteins
        fu = max(0.01, min(0.95, fu))
        
        # Base tissue partition coefficients
        kp_liver = 1.5 + logp * 0.5
        kp_brain = 0.2 + logp * 0.4  # BBB penetration
        kp_muscle = 0.5 + logp * 0.2
        kp_fat = 2.0 + logp * 1.0
        kp_lung = 1.0 + logp * 0.3
        kp_heart = 0.8 + logp * 0.2
        kp_kidney = 1.2 + logp * 0.3
        kp_gut = 1.0 + logp * 0.4
        kp_pancreas = 0.9 + logp * 0.25
        
        # Apply disease-specific tissue distribution adjustments
        if 'kp_brain' in self.pbpk_adjustments:
            kp_brain *= self.pbpk_adjustments['kp_brain']
        if 'bbb_factor' in self.pbpk_adjustments:
            kp_brain *= self.pbpk_adjustments['bbb_factor']
        if 'kp_liver' in self.pbpk_adjustments:
            kp_liver *= self.pbpk_adjustments['kp_liver']
        if 'kp_lung' in self.pbpk_adjustments:
            kp_lung *= self.pbpk_adjustments['kp_lung']
        if 'kp_heart' in self.pbpk_adjustments:
            kp_heart *= self.pbpk_adjustments['kp_heart']
        if 'kp_kidney' in self.pbpk_adjustments:
            kp_kidney *= self.pbpk_adjustments['kp_kidney']
        if 'kp_gut' in self.pbpk_adjustments:
            kp_gut *= self.pbpk_adjustments['kp_gut']
        if 'kp_pancreas' in self.pbpk_adjustments:
            kp_pancreas *= self.pbpk_adjustments['kp_pancreas']
        if 'tumor_penetration' in self.pbpk_adjustments:
            # Increase distribution to target tissue for oncology
            kp_liver *= self.pbpk_adjustments['tumor_penetration']
        if 'respiratory_distribution' in self.pbpk_adjustments:
            kp_lung *= self.pbpk_adjustments['respiratory_distribution']
        if 'immune_distribution' in self.pbpk_adjustments:
            # For immune-targeted diseases
            kp_liver *= self.pbpk_adjustments['immune_distribution']
        
        # Adjust based on binding affinity if available
        if binding_affinity and binding_affinity < -8.0:
            # Strong binding suggests high target affinity
            kp_liver *= 1.3  # Accumulates more in liver
            kp_brain *= 1.2  # Better brain penetration
        
        return {
            'ka': round(ka, 3),  # Absorption rate constant (h^-1)
            'vd': round(vd, 3),  # Volume of distribution (L/kg)
            'cl': round(cl, 3),  # Clearance (L/h/kg)
            'f': round(f, 3),  # Bioavailability
            'fu': round(fu, 3),  # Fraction unbound
            'kp_liver': round(kp_liver, 2),
            'kp_brain': round(kp_brain, 2),
            'kp_muscle': round(kp_muscle, 2),
            'kp_fat': round(kp_fat, 2),
            'kp_lung': round(kp_lung, 2),
            'kp_heart': round(kp_heart, 2),
            'kp_kidney': round(kp_kidney, 2),
            'kp_gut': round(kp_gut, 2),
            'kp_pancreas': round(kp_pancreas, 2),
            't_half': round(0.693 * vd / cl, 2),  # Half-life (hours)
            'disease_name': self.disease_name,
            'target_tissues': self.target_tissues
        }
    
    def _simulate_oral_dosing(
        self,
        dose_mg: float,
        duration_hours: float,
        adme: Dict
    ) -> Tuple[List[float], Dict[str, List[float]]]:
        """Simulate oral dosing with simplified exponential model"""
        
        # Time points (every 0.1 hours)
        num_points = int(duration_hours * 10)
        time_points = [i * duration_hours / num_points for i in range(num_points)]
        
        # Simplified oral absorption model (first-order kinetics)
        ka = adme['ka']
        ke = adme['cl'] / adme['vd']  # Elimination rate constant
        vd = adme['vd'] * self.human_params['body_weight']
        f = adme['f']
        
        # Analytical solution for oral one-compartment model
        plasma_conc = []
        liver_conc = []
        brain_conc = []
        kidney_conc = []
        
        for t in time_points:
            if ka != ke:
                # Standard one-compartment oral model
                c_plasma = (f * dose_mg * ka / (vd * (ka - ke))) * (math.exp(-ke * t) - math.exp(-ka * t))
            else:
                # Special case when ka == ke
                c_plasma = (f * dose_mg * ka * t / vd) * math.exp(-ke * t)
            
            c_plasma = max(0, c_plasma)  # Ensure non-negative
            
            # Tissue distribution based on partition coefficients
            c_liver = c_plasma * adme['kp_liver'] * (1 - math.exp(-0.3 * t))
            c_brain = c_plasma * adme['kp_brain'] * (1 - math.exp(-0.2 * t))
            c_kidney = c_plasma * 0.8  # Kidney follows plasma closely
            
            # Convert to ng/mL
            plasma_conc.append(c_plasma * 1000)
            liver_conc.append(c_liver * 1000)
            brain_conc.append(c_brain * 1000)
            kidney_conc.append(c_kidney * 1000)
        
        concentrations = {
            'plasma': plasma_conc,
            'liver': liver_conc,
            'kidney': kidney_conc,
            'brain': brain_conc,
            'target': liver_conc  # Assume liver as target
        }
        
        return time_points, concentrations
    
    def _simulate_iv_dosing(
        self,
        dose_mg: float,
        duration_hours: float,
        adme: Dict
    ) -> Tuple[List[float], Dict[str, List[float]]]:
        """Simulate IV dosing (bolus injection)"""
        
        num_points = int(duration_hours * 10)
        time_points = [i * duration_hours / num_points for i in range(num_points)]
        
        vd = adme['vd'] * self.human_params['body_weight']
        cl = adme['cl'] * self.human_params['body_weight']
        ke = cl / vd  # Elimination rate constant
        
        # One-compartment model for IV
        plasma_conc = []
        liver_conc = []
        brain_conc = []
        kidney_conc = []
        
        for t in time_points:
            c_plasma = (dose_mg / vd) * math.exp(-ke * t)
            c_liver = c_plasma * adme['kp_liver'] * (1 - math.exp(-0.3 * t))
            c_brain = c_plasma * adme['kp_brain'] * (1 - math.exp(-0.2 * t))
            c_kidney = c_plasma * 0.8
            
            plasma_conc.append(c_plasma * 1000)
            liver_conc.append(c_liver * 1000)
            brain_conc.append(c_brain * 1000)
            kidney_conc.append(c_kidney * 1000)
        
        concentrations = {
            'plasma': plasma_conc,
            'liver': liver_conc,
            'brain': brain_conc,
            'kidney': kidney_conc,
            'target': liver_conc
        }
        
        return time_points, concentrations
    
    def _calculate_pk_metrics(
        self,
        time_points: List[float],
        concentrations: Dict[str, List[float]],
        adme: Dict
    ) -> Dict:
        """Calculate pharmacokinetic metrics"""
        
        plasma_conc = concentrations['plasma']
        
        # Cmax: Maximum concentration
        cmax = max(plasma_conc)
        tmax = time_points[plasma_conc.index(cmax)]
        
        # AUC: Area under curve (trapezoidal rule)
        auc = 0
        for i in range(len(time_points) - 1):
            dt = time_points[i+1] - time_points[i]
            avg_conc = (plasma_conc[i] + plasma_conc[i+1]) / 2
            auc += avg_conc * dt
        
        # Time above threshold (assume threshold = 100 ng/mL for efficacy)
        dt = time_points[1] - time_points[0] if len(time_points) > 1 else 0
        time_above_threshold = sum(1 for c in plasma_conc if c > 100) * dt
        
        return {
            'cmax_ng_ml': round(cmax, 2),
            'tmax_hours': round(tmax, 2),
            'auc_ng_h_ml': round(auc, 2),
            't_half_hours': adme['t_half'],
            'time_above_threshold_hours': round(time_above_threshold, 2),
            'vd_l': round(adme['vd'] * self.human_params['body_weight'], 2),
            'clearance_l_h': round(adme['cl'] * self.human_params['body_weight'], 2)
        }
    
    def _assess_safety(self, pk_metrics: Dict, binding_affinity: Optional[float]) -> Dict:
        """Assess safety based on PK metrics"""
        
        cmax = pk_metrics['cmax_ng_ml']
        auc = pk_metrics['auc_ng_h_ml']
        
        # Safety thresholds (example values)
        safety_margin = "Good"
        warnings = []
        
        if cmax > 10000:
            safety_margin = "Caution"
            warnings.append("High peak concentration may cause acute toxicity")
        
        if auc > 100000:
            safety_margin = "Caution"
            warnings.append("High exposure may lead to accumulation")
        
        if binding_affinity and binding_affinity < -12.0:
            warnings.append("Very strong binding may cause off-target effects")
        
        if not warnings:
            warnings.append("No significant safety concerns identified")
        
        return {
            'safety_margin': safety_margin,
            'warnings': warnings,
            'therapeutic_window': "Within acceptable range" if safety_margin == "Good" else "Requires monitoring"
        }
    
    def analyze_repurposing_feasibility(
        self,
        drug_name: str,
        molecular_weight: float,
        logp: float,
        binding_affinity_kcal_mol: float,
        target_organ: str,
        dose_mg: float = 100.0,
        route: str = 'oral',
        original_indication: Optional[str] = None,
        new_indication: Optional[str] = None
    ) -> Dict:
        """
        Analyze drug repurposing feasibility based on PBPK and target binding
        
        Args:
            drug_name: Name of the drug
            molecular_weight: Molecular weight (g/mol)
            logp: Lipophilicity
            binding_affinity_kcal_mol: Docking binding affinity (kcal/mol, negative values)
            target_organ: Target organ/tissue (brain, liver, kidney, muscle)
            dose_mg: Proposed dose
            route: Administration route
            original_indication: Original disease indication
            new_indication: Proposed new indication
        
        Returns:
            Comprehensive feasibility assessment
        """
        logger.info(f"Analyzing repurposing feasibility for {drug_name} targeting {target_organ}")
        
        # Run PBPK simulation
        pbpk_result = self.simulate_drug_exposure(
            drug_name=drug_name,
            molecular_weight=molecular_weight,
            logp=logp,
            dose_mg=dose_mg,
            route=route,
            duration_hours=24.0,
            binding_affinity=binding_affinity_kcal_mol
        )
        
        # Convert binding affinity to Ki (nM) for concentration comparison
        # ΔG = -RT ln(Ki), therefore Ki = exp(ΔG/(-RT))
        # where ΔG in kcal/mol (negative for favorable binding), R = 0.001987 kcal/(mol·K), T = 298 K
        import math
        RT = 0.001987 * 298  # 0.592 kcal/mol at 298K
        # For binding affinity (negative values), Ki = exp(ΔG/RT) with NEGATIVE binding energy
        # This gives SMALL Ki for STRONG (negative) binding affinity
        Ki_M = math.exp(binding_affinity_kcal_mol / RT)  # Ki in Molar - corrected sign!
        Ki_nM = Ki_M * 1e9  # Convert to nM
        therapeutic_threshold_ng_ml = Ki_nM * (molecular_weight / 1000)  # Convert to ng/mL
        
        # Get target tissue concentration
        adme = pbpk_result['adme_parameters']
        pk_metrics = pbpk_result['pk_metrics']
        
        target_concentrations = {
            'brain': pbpk_result['brain_concentration_ng_ml'],
            'liver': pbpk_result['liver_concentration_ng_ml'],
            'plasma': pbpk_result['plasma_concentration_ng_ml'],
            'kidney': pbpk_result['plasma_concentration_ng_ml'],  # Approximation
            'muscle': [c * 0.5 for c in pbpk_result['plasma_concentration_ng_ml']]  # Approximation
        }
        
        target_conc_list = target_concentrations.get(target_organ.lower(), pbpk_result['plasma_concentration_ng_ml'])
        max_target_conc = max(target_conc_list) if target_conc_list else 0
        avg_target_conc = sum(target_conc_list) / len(target_conc_list) if target_conc_list else 0
        
        # Calculate exposure adequacy ratio
        exposure_ratio = max_target_conc / therapeutic_threshold_ng_ml if therapeutic_threshold_ng_ml > 0 else 0
        
        # Assess BBB penetration for brain targets
        bbb_penetration = None
        if target_organ.lower() == 'brain':
            brain_plasma_ratio = adme['kp_brain']
            if brain_plasma_ratio > 1.0:
                bbb_status = "Excellent"
                bbb_score = 100
            elif brain_plasma_ratio > 0.5:
                bbb_status = "Good"
                bbb_score = 75
            elif brain_plasma_ratio > 0.2:
                bbb_status = "Moderate"
                bbb_score = 50
            else:
                bbb_status = "Poor"
                bbb_score = 25
            
            bbb_penetration = {
                'brain_plasma_ratio': round(brain_plasma_ratio, 3),
                'status': bbb_status,
                'score': bbb_score,
                'message': f"Brain-to-plasma ratio: {brain_plasma_ratio:.2f}. LogP: {logp:.2f} (optimal BBB range: 1.5-2.5)"
            }
        
        # Determine feasibility verdict
        feasibility_score = 0
        feasibility_factors = []
        
        # Factor 1: Target reachability (40%)
        if target_organ.lower() == 'brain':
            if bbb_penetration:
                target_score = bbb_penetration['score']
                feasibility_factors.append(f"BBB Penetration: {bbb_penetration['status']} ({target_score}%)")
        else:
            # For non-brain targets, use available partition coefficient or estimate from concentration achieved
            target_kp = adme.get(f'kp_{target_organ.lower()}', None)
            if target_kp is None:
                # Estimate reachability from actual concentration ratio achieved
                if max_target_conc > 0 and pk_metrics['cmax_ng_ml'] > 0:
                    target_kp = max_target_conc / pk_metrics['cmax_ng_ml']
                else:
                    target_kp = 0.8  # Default: assume reasonable tissue distribution
            
            # Score based on partition coefficient
            if target_kp >= 1.5:
                target_score = 100
                reach_status = "Excellent tissue accumulation"
            elif target_kp >= 1.0:
                target_score = 85
                reach_status = "Good tissue distribution"
            elif target_kp >= 0.5:
                target_score = 70
                reach_status = "Moderate tissue distribution"
            elif target_kp >= 0.2:
                target_score = 50
                reach_status = "Fair tissue distribution"
            else:
                target_score = 30
                reach_status = "Limited tissue penetration"
            
            feasibility_factors.append(f"Tissue Distribution: {reach_status}")
        
        feasibility_score += target_score * 0.4
        
        # Factor 2: Concentration adequacy (40%)
        if exposure_ratio >= 10:
            conc_score = 100
            conc_status = "Excellent - 10x therapeutic threshold"
        elif exposure_ratio >= 5:
            conc_score = 85
            conc_status = "Very Good - 5x therapeutic threshold"
        elif exposure_ratio >= 2:
            conc_score = 70
            conc_status = "Good - 2x therapeutic threshold"
        elif exposure_ratio >= 1:
            conc_score = 50
            conc_status = "Adequate - meets therapeutic threshold"
        elif exposure_ratio >= 0.5:
            conc_score = 30
            conc_status = "Marginal - below therapeutic threshold"
        else:
            conc_score = 10
            conc_status = "Insufficient - far below therapeutic threshold"
        
        feasibility_factors.append(f"Target Concentration: {conc_status}")
        feasibility_score += conc_score * 0.4
        
        # Factor 3: Safety margin (20%)
        safety_ratio = pk_metrics['cmax_ng_ml'] / therapeutic_threshold_ng_ml if therapeutic_threshold_ng_ml > 0 else 0
        if safety_ratio < 50:
            safety_score = 100
            safety_status = "Excellent safety window"
        elif safety_ratio < 100:
            safety_score = 75
            safety_status = "Good safety window"
        elif safety_ratio < 200:
            safety_score = 50
            safety_status = "Moderate safety window"
        else:
            safety_score = 25
            safety_status = "Narrow safety window"
        
        feasibility_factors.append(f"Safety: {safety_status}")
        feasibility_score += safety_score * 0.2
        
        # Overall verdict
        if feasibility_score >= 80:
            verdict = "Highly Feasible"
            verdict_color = "success"
            recommendation = f"{drug_name} shows excellent potential for {new_indication or 'new indication'}. Target tissue exposure is adequate with good safety margins."
        elif feasibility_score >= 60:
            verdict = "Feasible with Optimization"
            verdict_color = "info"
            recommendation = f"{drug_name} is promising for {new_indication or 'new indication'}, but may benefit from dose optimization or formulation changes."
        elif feasibility_score >= 40:
            verdict = "Marginal Feasibility"
            verdict_color = "warning"
            recommendation = f"{drug_name} shows limited potential for {new_indication or 'new indication'}. Significant optimization or combination therapy may be required."
        else:
            verdict = "Low Feasibility"
            verdict_color = "error"
            recommendation = f"{drug_name} is unlikely to be effective for {new_indication or 'new indication'} at standard doses. Consider alternative candidates."
        
        # Dose adjustment recommendations
        dose_recommendations = []
        if exposure_ratio < 1.0 and exposure_ratio > 0:
            suggested_dose = dose_mg * (1.5 / exposure_ratio)
            suggested_dose = min(suggested_dose, dose_mg * 5)  # Cap at 5x increase
            dose_recommendations.append(f"Consider increasing dose to approximately {suggested_dose:.0f} mg to achieve therapeutic exposure")
        
        if exposure_ratio > 10:
            suggested_dose = dose_mg / 2
            dose_recommendations.append(f"Current dose may be excessive. Consider reducing to {suggested_dose:.0f} mg")
        
        if target_organ.lower() == 'brain' and bbb_penetration and bbb_penetration['score'] < 50:
            dose_recommendations.append("Poor BBB penetration. Consider: (1) Lipophilic prodrug, (2) Nanoparticle formulation, (3) Alternative route (intranasal)")
        
        if not dose_recommendations:
            dose_recommendations.append(f"Current dose of {dose_mg} mg appears appropriate")
        
        return {
            'success': True,
            'drug_name': drug_name,
            'original_indication': original_indication or 'Unknown',
            'new_indication': new_indication or 'Unknown',
            'target_organ': target_organ,
            'binding_affinity_kcal_mol': round(binding_affinity_kcal_mol, 2),
            'ki_nM': round(Ki_nM, 2),
            'therapeutic_threshold_ng_ml': round(therapeutic_threshold_ng_ml, 2),
            'max_target_concentration_ng_ml': round(max_target_conc, 2),
            'avg_target_concentration_ng_ml': round(avg_target_conc, 2),
            'exposure_adequacy_ratio': round(exposure_ratio, 2),
            'bbb_penetration': bbb_penetration,
            'feasibility_score': round(feasibility_score, 1),
            'feasibility_verdict': verdict,
            'verdict_color': verdict_color,
            'feasibility_factors': feasibility_factors,
            'recommendation': recommendation,
            'dose_mg': dose_mg,
            'route': route,
            'dose_recommendations': dose_recommendations,
            'pk_metrics': pk_metrics,
            'adme_parameters': adme,
            'safety_assessment': pbpk_result['safety_assessment'],
            'time_hours': pbpk_result['time_hours'],
            'plasma_concentration_ng_ml': pbpk_result['plasma_concentration_ng_ml'],
            'target_concentration_ng_ml': target_conc_list
        }


# Global instance
pbpk_simulator = PBPKSimulator()
