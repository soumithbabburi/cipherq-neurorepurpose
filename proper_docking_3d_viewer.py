"""
PROPER Docking 3D Viewer
Uses ACTUAL docking coordinates from AutoDock Vina/NVIDIA
Shows drug EXACTLY where it binds in the protein pocket
"""
import streamlit as st
import streamlit.components.v1 as components
import logging

logger = logging.getLogger(__name__)

def render_docking_3d_viewer(
    protein_pdb_content: str,
    ligand_sdf_content: str,
    drug_name: str,
    target_protein: str,
    binding_affinity: float,
    pose_number: int = 1,
    height: int = 700
):
    """
    Render 3D viewer showing drug EXACTLY where it docks in protein.
    Uses coordinates from docking SDF file - NOT separate loading!
    """
    
    # Escape strings for JavaScript
    protein_pdb = protein_pdb_content.replace('`', '\\`').replace('${', '\\${')
    ligand_sdf = ligand_sdf_content.replace('`', '\\`').replace('${', '\\${')
    
    html_code = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        <style>
            body {{ margin: 0; padding: 0; }}
            #container {{ width: 100%; height: {height}px; position: relative; }}
            .info-overlay {{
                position: absolute;
                top: 15px;
                left: 15px;
                background: rgba(255, 255, 255, 0.95);
                padding: 15px 20px;
                border-radius: 10px;
                box-shadow: 0 4px 12px rgba(0,0,0,0.15);
                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
                z-index: 1000;
            }}
            .info-title {{ 
                font-size: 16px; 
                font-weight: 600; 
                color: #1e293b; 
                margin-bottom: 8px;
            }}
            .info-item {{ 
                font-size: 13px; 
                color: #64748b; 
                margin: 3px 0;
            }}
            .affinity {{ 
                font-size: 20px; 
                font-weight: 700; 
                color: #059669; 
                margin-top: 8px;
            }}
            .controls {{
                position: absolute;
                bottom: 15px;
                left: 15px;
                background: rgba(255, 255, 255, 0.9);
                padding: 10px 15px;
                border-radius: 8px;
                font-size: 12px;
                color: #64748b;
            }}
        </style>
    </head>
    <body>
        <div id="container">
            <div class="info-overlay">
                <div class="info-title">🧬 Molecular Docking Complex</div>
                <div class="info-item"><strong>Drug:</strong> {drug_name}</div>
                <div class="info-item"><strong>Target:</strong> {target_protein}</div>
                <div class="info-item"><strong>Pose:</strong> #{pose_number}</div>
                <div class="affinity">{binding_affinity:.1f} kcal/mol</div>
            </div>
            <div class="controls">
                🖱️ Rotate: Left-drag | Zoom: Scroll | Pan: Right-drag | Reset: Double-click
            </div>
        </div>
        
        <script>
            let viewer = $3Dmol.createViewer(document.getElementById('container'), {{
                backgroundColor: 'white'
            }});
            
            // STEP 1: Add protein structure
            let proteinData = `{protein_pdb}`;
            viewer.addModel(proteinData, "pdb");
            
            // Style protein as cartoon (blue/gray)
            viewer.setStyle({{model: 0}}, {{
                cartoon: {{color: 'lightblue', opacity: 0.7}},
            }});
            
            // STEP 2: Add ligand from DOCKING SDF (has correct coordinates!)
            let ligandData = `{ligand_sdf}`;
            viewer.addModel(ligandData, "sdf");
            
            // Style ligand as ball-and-stick (colorful)
            viewer.setStyle({{model: 1}}, {{
                stick: {{
                    colorscheme: 'Jmol',
                    radius: 0.25
                }},
                sphere: {{
                    scale: 0.35,
                    colorscheme: 'Jmol'
                }}
            }});
            
            // Add surface to protein for better visualization
            viewer.addSurface($3Dmol.SurfaceType.VDW, {{
                opacity: 0.15,
                color: 'lightgray'
            }}, {{model: 0}});
            
            // CRITICAL: Center on LIGAND (drug)
            // This ensures we see the binding pocket, not the whole protein
            viewer.zoomTo({{model: 1}});
            
            // Zoom out slightly to see protein context around binding site
            viewer.zoom(0.7, 1000);
            
            // Render
            viewer.render();
            
            // Auto-rotate on double-click
            let rotating = false;
            document.getElementById('container').addEventListener('dblclick', function() {{
                rotating = !rotating;
                if (rotating) {{
                    viewer.spin(true);
                }} else {{
                    viewer.spin(false);
                }}
            }});
            
            console.log('3D viewer loaded: {drug_name} → {target_protein}');
            console.log('Ligand centered in binding pocket');
        </script>
    </body>
    </html>
    """
    
    components.html(html_code, height=height)
    
    st.caption(f"💡 **The drug is positioned at the exact docking coordinates.** Double-click to rotate.")
    
    logger.info(f"✅ Rendered docking complex: {drug_name} in {target_protein} binding pocket")
    return True


# Alternative: Side-by-side comparison
def render_docking_comparison(
    protein_pdb: str,
    original_sdf: str,
    optimized_sdf: str,
    drug_name: str,
    target_protein: str,
    original_affinity: float,
    optimized_affinity: float,
    height: int = 600
):
    """
    Show side-by-side: Original drug vs Optimized drug in same protein pocket.
    Both use actual docking coordinates.
    """
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown(f"### Original: {drug_name}")
        st.markdown(f"**Affinity:** {original_affinity:.1f} kcal/mol")
        render_docking_3d_viewer(
            protein_pdb, original_sdf, drug_name, target_protein, 
            original_affinity, 1, height
        )
    
    with col2:
        st.markdown(f"### Optimized: {drug_name}")
        st.markdown(f"**Affinity:** {optimized_affinity:.1f} kcal/mol")
        st.markdown(f"**Change:** {optimized_affinity - original_affinity:+.1f} kcal/mol")
        render_docking_3d_viewer(
            protein_pdb, optimized_sdf, f"{drug_name} (optimized)", target_protein,
            optimized_affinity, 1, height
        )
    
    return True
