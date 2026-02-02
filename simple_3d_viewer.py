"""
Simple 3D Viewer Stub Module
Fallback implementation when full 3D viewer is not available
"""

import streamlit as st

def create_simple_3d_viewer(molecule_data=None, protein_data=None, title="3D Molecular Viewer"):
    """
    Fallback 3D viewer using basic visualization
    
    Args:
        molecule_data: Drug molecule structure data
        protein_data: Protein structure data
        title: Viewer title
    """
    st.info(f"🔬 {title}")
    
    if molecule_data:
        st.write("**Drug Structure:**")
        if isinstance(molecule_data, dict):
            st.json(molecule_data)
        else:
            st.text(str(molecule_data))
    
    if protein_data:
        st.write("**Protein Structure:**")
        if isinstance(protein_data, dict):
            st.json(protein_data)
        else:
            st.text(str(protein_data))
    
    st.warning("⚠️ Advanced 3D visualization requires additional dependencies. Showing simplified view.")
    
    return None

def render_3d_structure(structure_data, viewer_type="molecule"):
    """Fallback renderer"""
    st.write(f"**{viewer_type.title()} Structure Data:**")
    st.json(structure_data if isinstance(structure_data, dict) else {"data": str(structure_data)})
    return None

# Ensure module can be imported
__all__ = ['create_simple_3d_viewer', 'render_3d_structure']
