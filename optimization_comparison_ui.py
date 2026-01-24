#!/usr/bin/env python3
"""
OPTIMIZATION COMPARISON UI
Displays before/after scores and real improvements
"""

import streamlit as st
import plotly.graph_objects as go
from typing import Dict, List
from real_molecular_optimizer import OptimizationResult


def display_score_card(title: str, score: float, subtitle: str = "", delta: float = None, tooltip: str = ""):
    """Display a score card with optional improvement delta - CLEAN VERSION"""
    
    # Determine color based on score
    if score >= 80:
        score_text = f":green[{score:.1f}%]"
    elif score >= 60:
        score_text = f":orange[{score:.1f}%]"
    else:
        score_text = f":red[{score:.1f}%]"
    
    # Display using Streamlit metrics
    if delta is not None:
        st.metric(
            label=title,
            value=f"{score:.1f}%",
            delta=f"{delta:+.1f}%",
            help=tooltip if tooltip else None
        )
        if subtitle:
            st.caption(subtitle)
    else:
        st.metric(
            label=title,
            value=f"{score:.1f}%",
            help=tooltip if tooltip else None
        )
        if subtitle:
            st.caption(subtitle)


def display_property_comparison(original_props: Dict, optimized_props: Dict, changes: Dict):
    """Display side-by-side property comparison - CLEAN VERSION"""
    
    # Key properties to show
    key_props = [
        ('LogP', 'Lipophilicity', ''),
        ('TPSA', 'Polar Surface Area', 'A^2'),
        ('MW', 'Molecular Weight', 'g/mol'),
        ('BBB_Score', 'BBB Penetration', '%'),
        ('CNS_MPO', 'CNS MPO Score', '/6'),
        ('DrugLikeness', 'Drug-Likeness', '%'),
    ]
    
    st.markdown("### Property Comparison")
    
    for prop_key, prop_name, unit in key_props:
        if prop_key in original_props and prop_key in optimized_props:
            orig_val = original_props[prop_key]
            opt_val = optimized_props[prop_key]
            change = changes.get(prop_key, 0)
            
            # Create three-column layout
            col1, col2, col3, col4 = st.columns([2, 1.5, 1.5, 1.5])
            
            with col1:
                st.markdown(f"**{prop_name}**")
            
            with col2:
                st.text("Original")
                st.markdown(f"`{orig_val:.2f} {unit}`")
            
            with col3:
                st.text("Optimized")
                st.markdown(f"`{opt_val:.2f} {unit}`")
            
            with col4:
                # Determine if change is positive
                if prop_key == 'TPSA':
                    better = change < 0
                elif prop_key in ['LogP', 'BBB_Score', 'CNS_MPO', 'DrugLikeness']:
                    better = change > 0
                else:
                    better = abs(change) < 50
                
                change_symbol = "+" if change > 0 else ""
                if better:
                    st.markdown(f":green[{change_symbol}{change:.2f} {unit}]")
                else:
                    st.markdown(f":red[{change_symbol}{change:.2f} {unit}]")
            
            st.divider()


def display_optimization_results(results: List[OptimizationResult], drug_name: str):
    """Display complete optimization results with before/after comparison"""
    
    if not results or not results[0].success:
        st.error("No successful optimizations found")
        if results and results[0].error_message:
            st.error(f"Error: {results[0].error_message}")
        return
    
    # Header
    st.markdown(f"## Optimization Results: {drug_name}")
    st.markdown("---")
    
    # Show original score first
    st.markdown("### Initial Analysis")
    original_score = results[0].original_score
    original_props = results[0].original_properties
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        display_score_card(
            "Overall Score",
            original_score,
            "Before optimization"
        )
    
    with col2:
        bbb_score = original_props.get('BBB_Score', 0)
        display_score_card(
            "BBB Penetration",
            bbb_score,
            "Blood-brain barrier"
        )
    
    with col3:
        cns_mpo = original_props.get('CNS_MPO', 0) / 6.0 * 100
        display_score_card(
            "CNS Properties",
            cns_mpo,
            "CNS MPO normalized"
        )
    
    st.markdown("---")
    
    # Show best optimization
    best_result = results[0]
    
    st.markdown(f"### Best Optimization: {best_result.modification_type}")
    
    # Create optimization explanation tooltip
    mod_type = best_result.modification_type
    if "N-methylation" in mod_type:
        optimization_method = "Added methyl group to nitrogen atoms to increase lipophilicity and BBB penetration"
    elif "Ethylation" in mod_type:
        optimization_method = "Added ethyl groups to increase lipophilicity while maintaining hydrogen bonding capacity"
    elif "Polarity reduction" in mod_type or "dehydroxylation" in mod_type:
        optimization_method = "Removed hydroxyl groups to reduce polarity and improve BBB penetration"
    elif "Fluorination" in mod_type:
        optimization_method = "Added fluorine atoms to modulate metabolism and improve lipophilicity"
    else:
        optimization_method = "Chemical modification to optimize CNS drug properties"
    
    # Before/After Scores
    col1, col2, col3 = st.columns(3)
    
    with col1:
        display_score_card(
            "Optimized Score",
            best_result.optimized_score,
            "After optimization",
            delta=best_result.score_improvement,
            tooltip=f"Overall improvement via {optimization_method}"
        )
    
    with col2:
        opt_bbb = best_result.optimized_properties.get('BBB_Score', 0)
        orig_bbb = best_result.original_properties.get('BBB_Score', 0)
        display_score_card(
            "BBB Penetration",
            opt_bbb,
            "Optimized",
            delta=opt_bbb - orig_bbb,
            tooltip=f"BBB improvement via {optimization_method}"
        )
    
    with col3:
        opt_cns = best_result.optimized_properties.get('CNS_MPO', 0) / 6.0 * 100
        orig_cns = best_result.original_properties.get('CNS_MPO', 0) / 6.0 * 100
        display_score_card(
            "CNS Properties",
            opt_cns,
            "Optimized",
            delta=opt_cns - orig_cns,
            tooltip=f"CNS MPO score enhanced via {optimization_method}"
        )
    
    # Confidence boost message
    st.info(f"**{best_result.confidence_boost}**")
    
    # Detailed property comparison
    st.markdown("---")
    display_property_comparison(
        best_result.original_properties,
        best_result.optimized_properties,
        best_result.property_changes
    )
    
    # Chemical structures
    st.markdown("---")
    st.markdown("### Chemical Structures")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Original Structure**")
        st.code(best_result.original_smiles, language="text")
    
    with col2:
        st.markdown("**Optimized Structure**")
        st.code(best_result.optimized_smiles, language="text")
    
    # Show all other optimization attempts
    if len(results) > 1:
        with st.expander(f"View All {len(results)} Optimization Strategies"):
            for i, result in enumerate(results, 1):
                st.markdown(f"**{i}. {result.modification_type}**")
                st.write(f"Score: {result.original_score:.1f}% → {result.optimized_score:.1f}% ({result.score_improvement:+.1f}%)")
                st.write(f"SMILES: `{result.optimized_smiles}`")
                st.markdown("---")


def create_score_comparison_chart(results: List[OptimizationResult]) -> go.Figure:
    """Create visual comparison chart of all optimization attempts"""
    
    if not results:
        return None
    
    modifications = [r.modification_type for r in results]
    original_scores = [r.original_score for r in results]
    optimized_scores = [r.optimized_score for r in results]
    
    fig = go.Figure()
    
    # Original scores
    fig.add_trace(go.Bar(
        name='Original',
        x=modifications,
        y=original_scores,
        marker_color='#6c757d'
    ))
    
    # Optimized scores
    fig.add_trace(go.Bar(
        name='Optimized',
        x=modifications,
        y=optimized_scores,
        marker_color='#28a745'
    ))
    
    fig.update_layout(
        title='Optimization Strategy Comparison',
        xaxis_title='Modification Type',
        yaxis_title='Score (%)',
        barmode='group',
        height=400,
        yaxis=dict(range=[0, 100])
    )
    
    return fig
