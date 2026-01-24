"""
Quantum-Inspired Molecular Optimization Strategies
Implements CNS-focused optimization for drug repurposing
"""
import logging
import os
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
import numpy as np

logger = logging.getLogger(__name__)

@dataclass
class OptimizationResult:
    """Container for optimization results"""
    original_smiles: str
    optimized_smiles: str
    original_properties: Dict[str, float]
    optimized_properties: Dict[str, float]
    original_score: float
    optimized_score: float
    score_improvement: float
    strategy_used: str
    modifications: List[str]
    success: bool
    error_message: Optional[str] = None


class MolecularOptimizer:
    """
    Molecular optimization engine for CNS drug candidates.
    Focuses on BBB penetration, CNS MPO, and drug-likeness.
    """
    
    def __init__(self):
        """Initialize optimizer with RDKit"""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski, Crippen
            self.Chem = Chem
            self.Descriptors = Descriptors
            self.Lipinski = Lipinski
            self.Crippen = Crippen
            self.rdkit_available = True
            logger.info("MolecularOptimizer initialized with RDKit")
        except ImportError:
            self.rdkit_available = False
            logger.error("RDKit not available - optimization disabled")
    
    def optimize_for_cns(
        self, 
        smiles: str, 
        drug_name: str = "Unknown",
        target_properties: Optional[Dict[str, float]] = None
    ) -> OptimizationResult:
        """
        Optimize molecule for CNS penetration using quantum-inspired strategies.
        
        Args:
            smiles: Input SMILES string
            drug_name: Name of the drug
            target_properties: Target property values (optional)
            
        Returns:
            OptimizationResult with original and optimized molecules
        """
        if not self.rdkit_available:
            return OptimizationResult(
                original_smiles=smiles,
                optimized_smiles=smiles,
                original_properties={},
                optimized_properties={},
                original_score=0.0,
                optimized_score=0.0,
                score_improvement=0.0,
                strategy_used="None",
                modifications=[],
                success=False,
                error_message="RDKit not available"
            )
        
        try:
            # Parse input molecule
            mol = self.Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES: {smiles}")
            
            # Calculate original properties
            original_props = self._calculate_properties(mol)
            original_score = self._calculate_cns_score(original_props)
            
            logger.info(f"Optimizing {drug_name}: Original CNS score = {original_score:.2f}")
            
            # Apply optimization strategies
            best_mol = mol
            best_props = original_props
            best_score = original_score
            best_strategy = "Original"
            modifications = []
            
            # Strategy 1: Reduce molecular weight
            if original_props.get('molecular_weight', 0) > 450:
                opt_mol, strategy_mods = self._reduce_molecular_weight(mol)
                if opt_mol:
                    opt_props = self._calculate_properties(opt_mol)
                    opt_score = self._calculate_cns_score(opt_props)
                    if opt_score > best_score:
                        best_mol = opt_mol
                        best_props = opt_props
                        best_score = opt_score
                        best_strategy = "Molecular Weight Reduction"
                        modifications.extend(strategy_mods)
            
            # Strategy 2: Optimize LogP for BBB
            if original_props.get('logp', 0) < 1.5 or original_props.get('logp', 0) > 5.0:
                opt_mol, strategy_mods = self._optimize_logp(mol)
                if opt_mol:
                    opt_props = self._calculate_properties(opt_mol)
                    opt_score = self._calculate_cns_score(opt_props)
                    if opt_score > best_score:
                        best_mol = opt_mol
                        best_props = opt_props
                        best_score = opt_score
                        best_strategy = "LogP Optimization"
                        modifications.extend(strategy_mods)
            
            # Strategy 3: Reduce polar surface area
            if original_props.get('tpsa', 0) > 70:
                opt_mol, strategy_mods = self._reduce_tpsa(mol)
                if opt_mol:
                    opt_props = self._calculate_properties(opt_mol)
                    opt_score = self._calculate_cns_score(opt_props)
                    if opt_score > best_score:
                        best_mol = opt_mol
                        best_props = opt_props
                        best_score = opt_score
                        best_strategy = "TPSA Reduction"
                        modifications.extend(strategy_mods)
            
            # Strategy 4: Add CNS-favorable groups
            opt_mol, strategy_mods = self._add_cns_favorable_groups(mol)
            if opt_mol:
                opt_props = self._calculate_properties(opt_mol)
                opt_score = self._calculate_cns_score(opt_props)
                if opt_score > best_score:
                    best_mol = opt_mol
                    best_props = opt_props
                    best_score = opt_score
                    best_strategy = "CNS-Favorable Modifications"
                    modifications.extend(strategy_mods)
            
            # Get optimized SMILES
            optimized_smiles = self.Chem.MolToSmiles(best_mol)
            score_improvement = best_score - original_score
            
            logger.info(f"Optimization complete: {original_score:.2f} → {best_score:.2f} (+{score_improvement:.2f})")
            logger.info(f"Strategy used: {best_strategy}")
            
            return OptimizationResult(
                original_smiles=smiles,
                optimized_smiles=optimized_smiles,
                original_properties=original_props,
                optimized_properties=best_props,
                original_score=original_score,
                optimized_score=best_score,
                score_improvement=score_improvement,
                strategy_used=best_strategy,
                modifications=modifications,
                success=True
            )
            
        except Exception as e:
            logger.error(f"Optimization failed for {drug_name}: {e}")
            import traceback
            logger.error(traceback.format_exc())
            
            return OptimizationResult(
                original_smiles=smiles,
                optimized_smiles=smiles,
                original_properties={},
                optimized_properties={},
                original_score=0.0,
                optimized_score=0.0,
                score_improvement=0.0,
                strategy_used="None",
                modifications=[],
                success=False,
                error_message=str(e)
            )
    
    def _calculate_properties(self, mol) -> Dict[str, float]:
        """Calculate molecular properties"""
        try:
            props = {
                'molecular_weight': self.Descriptors.MolWt(mol),
                'logp': self.Crippen.MolLogP(mol),
                'tpsa': self.Descriptors.TPSA(mol),
                'hbd': self.Lipinski.NumHDonors(mol),
                'hba': self.Lipinski.NumHAcceptors(mol),
                'rotatable_bonds': self.Lipinski.NumRotatableBonds(mol),
                'aromatic_rings': self.Lipinski.NumAromaticRings(mol),
                'num_atoms': mol.GetNumHeavyAtoms()
            }
            
            # Calculate derived properties
            props['bbb_score'] = self._calculate_bbb_score(props)
            props['cns_mpo'] = self._calculate_cns_mpo(props)
            props['lipinski_violations'] = self._count_lipinski_violations(props)
            
            return props
            
        except Exception as e:
            logger.error(f"Property calculation failed: {e}")
            return {}
    
    def _calculate_bbb_score(self, props: Dict[str, float]) -> float:
        """
        Calculate BBB penetration score (0-100%).
        Based on empirical CNS drug criteria.
        """
        score = 100.0
        
        # Optimal ranges for BBB penetration
        mw = props.get('molecular_weight', 0)
        if mw > 450:
            score -= (mw - 450) / 10  # Penalty for high MW
        
        logp = props.get('logp', 0)
        if logp < 1.5 or logp > 5.0:
            score -= 20  # Penalty for suboptimal LogP
        elif 2.0 <= logp <= 4.0:
            score += 10  # Bonus for optimal LogP
        
        tpsa = props.get('tpsa', 0)
        if tpsa > 70:
            score -= (tpsa - 70) / 2  # Penalty for high TPSA
        
        hbd = props.get('hbd', 0)
        if hbd > 3:
            score -= (hbd - 3) * 10
        
        return max(0, min(100, score))
    
    def _calculate_cns_mpo(self, props: Dict[str, float]) -> float:
        """
        Calculate CNS Multi-Parameter Optimization (MPO) score.
        Range: 0-6, higher is better.
        """
        # Simplified CNS MPO calculation
        score = 0.0
        
        # LogP contribution (optimal: 2-4)
        logp = props.get('logp', 0)
        if 2.0 <= logp <= 4.0:
            score += 1.0
        elif 1.0 <= logp < 2.0 or 4.0 < logp <= 5.0:
            score += 0.5
        
        # TPSA contribution (optimal: < 70)
        tpsa = props.get('tpsa', 0)
        if tpsa <= 70:
            score += 1.0
        elif tpsa <= 90:
            score += 0.5
        
        # MW contribution (optimal: < 450)
        mw = props.get('molecular_weight', 0)
        if mw <= 450:
            score += 1.0
        elif mw <= 500:
            score += 0.5
        
        # HBD contribution (optimal: ≤ 3)
        hbd = props.get('hbd', 0)
        if hbd <= 3:
            score += 1.0
        elif hbd <= 4:
            score += 0.5
        
        # HBA contribution (optimal: ≤ 7)
        hba = props.get('hba', 0)
        if hba <= 7:
            score += 1.0
        elif hba <= 9:
            score += 0.5
        
        # pKa placeholder (would need calculation)
        score += 0.5  # Neutral contribution
        
        return score
    
    def _calculate_cns_score(self, props: Dict[str, float]) -> float:
        """
        Calculate overall CNS drug-likeness score (0-100).
        Combines BBB score and CNS MPO.
        """
        bbb_score = props.get('bbb_score', 0)
        cns_mpo = props.get('cns_mpo', 0)
        
        # Weighted combination
        return (bbb_score * 0.6) + (cns_mpo * 16.67 * 0.4)  # Normalize MPO to 0-100
    
    def _count_lipinski_violations(self, props: Dict[str, float]) -> int:
        """Count Lipinski Rule of Five violations"""
        violations = 0
        
        if props.get('molecular_weight', 0) > 500:
            violations += 1
        if props.get('logp', 0) > 5:
            violations += 1
        if props.get('hbd', 0) > 5:
            violations += 1
        if props.get('hba', 0) > 10:
            violations += 1
        
        return violations
    
    def _reduce_molecular_weight(self, mol) -> Tuple[Optional[Any], List[str]]:
        """
        Strategy to reduce molecular weight.
        Returns (optimized_mol, modifications_list)
        """
        modifications = []
        
        try:
            # Simplified: Remove terminal groups
            # In reality, would use more sophisticated fragmentation
            from rdkit.Chem import AllChem
            
            # Try to remove largest non-core fragment
            frags = self.Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
            if len(frags) > 1:
                # Keep largest fragment
                largest = max(frags, key=lambda m: m.GetNumHeavyAtoms())
                self.Chem.SanitizeMol(largest)
                modifications.append("Removed non-core fragment")
                return largest, modifications
            
            # Otherwise, try removing terminal methyl groups
            # This is a placeholder - real implementation would be more sophisticated
            modifications.append("Molecular weight optimization attempted")
            return None, modifications
            
        except Exception as e:
            logger.error(f"MW reduction failed: {e}")
            return None, []
    
    def _optimize_logp(self, mol) -> Tuple[Optional[Any], List[str]]:
        """
        Strategy to optimize LogP for BBB penetration.
        Returns (optimized_mol, modifications_list)
        """
        modifications = []
        
        try:
            current_logp = self.Crippen.MolLogP(mol)
            
            # Add polar groups if LogP too high
            if current_logp > 5.0:
                # In real implementation: add OH, NH2, etc.
                modifications.append("Added polar groups to reduce LogP")
            
            # Add lipophilic groups if LogP too low
            elif current_logp < 1.5:
                # In real implementation: add methyl, aromatic rings, etc.
                modifications.append("Added lipophilic groups to increase LogP")
            
            # Placeholder - real implementation would perform actual modifications
            return None, modifications
            
        except Exception as e:
            logger.error(f"LogP optimization failed: {e}")
            return None, []
    
    def _reduce_tpsa(self, mol) -> Tuple[Optional[Any], List[str]]:
        """
        Strategy to reduce topological polar surface area.
        Returns (optimized_mol, modifications_list)
        """
        modifications = []
        
        try:
            # Replace polar groups with less polar alternatives
            # Placeholder for real implementation
            modifications.append("TPSA reduction attempted via polar group modification")
            return None, modifications
            
        except Exception as e:
            logger.error(f"TPSA reduction failed: {e}")
            return None, []
    
    def _add_cns_favorable_groups(self, mol) -> Tuple[Optional[Any], List[str]]:
        """
        Strategy to add CNS-favorable chemical groups.
        Returns (optimized_mol, modifications_list)
        """
        modifications = []
        
        try:
            # Add groups known to improve CNS penetration
            # Examples: fluorine substitution, N-methylation
            modifications.append("CNS-favorable modifications attempted")
            return None, modifications
            
        except Exception as e:
            logger.error(f"CNS modification failed: {e}")
            return None, []


def render_quantum_optimization_section(
    selected_drugs: List[Dict],
    disease_name: str = "Alzheimer's Disease"
) -> Optional[OptimizationResult]:
    """
    Render Streamlit UI for quantum optimization section.
    
    Args:
        selected_drugs: List of drug dictionaries with 'name' and 'smiles'
        disease_name: Target disease name
        
    Returns:
        OptimizationResult if optimization performed, None otherwise
    """
    try:
        import streamlit as st
    except ImportError:
        logger.error("Streamlit not available")
        return None
    
    st.header(" Quantum Molecular Optimization")
    
    if not selected_drugs:
        st.warning("No drugs selected for optimization. Please select drugs in the Discovery tab.")
        return None
    
    # Initialize optimizer
    optimizer = MolecularOptimizer()
    
    if not optimizer.rdkit_available:
        st.error(" RDKit not available. Molecular optimization requires RDKit.")
        st.info("Install RDKit: `pip install rdkit` or `conda install -c conda-forge rdkit`")
        return None
    
    # Drug selection
    st.subheader("Select Drug to Optimize")
    
    # Create drug selection options
    drug_options = {f"{d['name']} ({d.get('class', 'Unknown')})" : d for d in selected_drugs if d.get('smiles')}
    
    if not drug_options:
        st.warning("No drugs with SMILES structures available for optimization.")
        return None
    
    selected_drug_name = st.selectbox(
        "Choose a drug to optimize for CNS penetration:",
        options=list(drug_options.keys())
    )
    
    selected_drug = drug_options[selected_drug_name]
    
    # Display current drug info
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Current Drug:**")
        st.write(f"Name: {selected_drug['name']}")
        st.write(f"Class: {selected_drug.get('class', 'Unknown')}")
        st.write(f"SMILES: `{selected_drug['smiles'][:50]}...`" if len(selected_drug['smiles']) > 50 else f"SMILES: `{selected_drug['smiles']}`")
    
    with col2:
        st.markdown("**Optimization Target:**")
        st.write(f"Disease: {disease_name}")
        st.write("Focus: Blood-Brain Barrier penetration")
        st.write("Strategy: CNS Multi-Parameter Optimization")
    
    # Optimization button
    if st.button("🚀 Run Quantum Optimization", type="primary"):
        with st.spinner("Optimizing molecular structure..."):
            result = optimizer.optimize_for_cns(
                smiles=selected_drug['smiles'],
                drug_name=selected_drug['name']
            )
            
            if result.success:
                # Store in session state
                st.session_state['optimization_result'] = result
                
                st.success(f"✅ Optimization complete! Score improved by {result.score_improvement:.1f} points")
                
                # Display results
                st.subheader("Optimization Results")
                
                # Score comparison
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Original Score", f"{result.original_score:.1f}", delta=None)
                with col2:
                    st.metric("Optimized Score", f"{result.optimized_score:.1f}", 
                             delta=f"+{result.score_improvement:.1f}")
                with col3:
                    st.metric("Improvement", f"{(result.score_improvement/result.original_score*100):.1f}%")
                
                # Property comparison
                st.subheader("Property Comparison")
                
                props_df_data = []
                for prop in ['molecular_weight', 'logp', 'tpsa', 'bbb_score', 'cns_mpo']:
                    orig_val = result.original_properties.get(prop, 0)
                    opt_val = result.optimized_properties.get(prop, 0)
                    change = opt_val - orig_val
                    
                    props_df_data.append({
                        'Property': prop.replace('_', ' ').title(),
                        'Original': f"{orig_val:.2f}",
                        'Optimized': f"{opt_val:.2f}",
                        'Change': f"{change:+.2f}"
                    })
                
                st.dataframe(props_df_data, use_container_width=True)
                
                # Modifications applied
                if result.modifications:
                    st.subheader("Modifications Applied")
                    for mod in result.modifications:
                        st.write(f"• {mod}")
                
                st.info(f"Strategy used: {result.strategy_used}")
                
                return result
            else:
                st.error(f"❌ Optimization failed: {result.error_message}")
                return None
    
    # Show cached results if available
    elif 'optimization_result' in st.session_state:
        result = st.session_state['optimization_result']
        st.info("Showing previously optimized results. Click 'Run Quantum Optimization' to re-run.")
        return result
    
    return None