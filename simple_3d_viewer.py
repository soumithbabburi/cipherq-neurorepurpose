"""
SIMPLE 3D VIEWER - Bulletproof py3Dmol implementation with zoom fix
"""

import streamlit as st
import streamlit.components.v1 as components
import json
import base64
import os
import glob
import logging

logger = logging.getLogger(__name__)

def create_simple_3d_viewer(drug_name: str, target_protein: str):
    """
    Create a bulletproof 3D viewer with proper zoom handling
    """
    
    try:
        # Get DiffDock poses
        diffdock_dir = f"./diffdock_output/{drug_name}"
        ligand_files = glob.glob(f"{diffdock_dir}/pose_*.sdf")[:3]  # Limit to 3 poses
        
        if not ligand_files:
            st.warning(f"No DiffDock poses found for {drug_name}")
            return False
            
        # Get protein PDB file - first check diffdock output, then pdb_cache
        protein_data = None
        
        # Priority 1: Check diffdock output folder for protein.pdb (saved during docking)
        protein_in_diffdock = f"{diffdock_dir}/protein.pdb"
        if os.path.exists(protein_in_diffdock):
            with open(protein_in_diffdock, 'r') as f:
                protein_data = f.read()
            logger.info(f"Using protein from docking output: {protein_in_diffdock}")
        
        # Priority 2: Check pdb_cache
        if not protein_data or len(protein_data) < 500:
            pdb_files = glob.glob(f"./pdb_cache/*{target_protein}*.pdb")
            if not pdb_files:
                pdb_files = glob.glob(f"./pdb_cache/*.pdb")
            
            if pdb_files:
                with open(pdb_files[0], 'r') as f:
                    protein_data = f.read()
                logger.info(f"Using protein from cache: {pdb_files[0]}")
        
        if not protein_data or len(protein_data) < 100:
            st.error(f"No protein structure found for {target_protein}")
            return False
        
        protein_b64 = base64.b64encode(protein_data.encode()).decode()
        
        # Read ligand poses and encode
        ligand_data = []
        for sdf_file in ligand_files:
            with open(sdf_file, 'r') as f:
                sdf_content = f.read()
            if sdf_content.strip():
                sdf_b64 = base64.b64encode(sdf_content.encode()).decode()
                ligand_data.append(sdf_b64)
        
        if not ligand_data:
            st.error("No valid ligand poses found")
            return False
        
        # Create viewer with stable zoom
        viewer_id = f"viewer_{drug_name}_{target_protein}".replace(' ', '_').replace(':', '_').replace('-', '_')
        
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/1.8.0/3Dmol-min.js"></script>
            <style>
                #{viewer_id} {{
                    width: 100% !important;
                    height: 500px !important;
                    background: white !important;
                    border: 2px solid #4CAF50 !important;
                    border-radius: 8px !important;
                    position: relative !important;
                }}
            </style>
        </head>
        <body>
            <div id="{viewer_id}" style="width: 100%; height: 500px; background: white;"></div>
            
            <script>
                console.log("[3D Viewer] Initializing...");
                
                var globalViewer = null;
                var modelsLoaded = false;
                
                function initViewer() {{
                    if (typeof $3Dmol === 'undefined') {{
                        setTimeout(initViewer, 100);
                        return;
                    }}
                    
                    try {{
                        console.log("[3D Viewer] Creating viewer instance...");
                        
                        // Create viewer with config
                        globalViewer = $3Dmol.createViewer('{viewer_id}', {{
                            backgroundColor: 'white',
                            antialias: true,
                            cartoonQuality: 10
                        }});
                        
                        if (!globalViewer) {{
                            throw new Error("Failed to create viewer");
                        }}
                        
                        // Load protein
                        console.log("[3D Viewer] Loading protein...");
                        var proteinData = atob('{protein_b64}');
                        globalViewer.addModel(proteinData, "pdb");
                        
                        // Load ligands
                        var ligands = {json.dumps(ligand_data)};
                        var colors = ['#FF0000', '#0000FF', '#00FF00'];
                        
                        console.log("[3D Viewer] Loading " + ligands.length + " ligand poses...");
                        for (var i = 0; i < ligands.length; i++) {{
                            var ligandData = atob(ligands[i]);
                            globalViewer.addModel(ligandData, "sdf");
                        }}
                        
                        console.log("[3D Viewer] Applying styles...");
                        
                        // Style protein (model 0)
                        globalViewer.setStyle({{model: 0}}, {{
                            cartoon: {{
                                color: 'lightblue',
                                opacity: 0.8,
                                thickness: 0.5
                            }}
                        }});
                        
                        // Style ligands (models 1, 2, 3)
                        for (var i = 0; i < ligands.length; i++) {{
                            globalViewer.setStyle({{model: i + 1}}, {{
                                stick: {{
                                    color: colors[i],
                                    radius: 0.3,
                                    opacity: 1.0
                                }}
                            }});
                        }}
                        
                        modelsLoaded = true;
                        
                        // Initial render and zoom
                        globalViewer.render();
                        globalViewer.zoomTo();
                        
                        // Zoom out slightly so everything stays visible
                        setTimeout(function() {{
                            if (globalViewer) {{
                                globalViewer.zoom(0.85, 1000);
                                globalViewer.render();
                                console.log("[3D Viewer] ✅ Viewer ready!");
                            }}
                        }}, 100);
                        
                        // Keep viewer responsive
                        setInterval(function() {{
                            if (globalViewer && modelsLoaded) {{
                                try {{
                                    globalViewer.render();
                                }} catch(e) {{
                                    // Ignore render errors
                                }}
                            }}
                        }}, 5000);
                        
                    }} catch(e) {{
                        console.error("[3D Viewer] ERROR:", e);
                        var elem = document.getElementById('{viewer_id}');
                        if (elem) {{
                            elem.innerHTML = '<div style="padding: 40px; text-align: center; color: #666;">' +
                                '⚠️ 3D visualization error<br/>' +
                                '<small>' + e.message + '</small><br/><br/>' +
                                '<button onclick="location.reload()" style="padding: 10px 20px; background: #4CAF50; color: white; border: none; border-radius: 4px; cursor: pointer;">Reload Page</button>' +
                                '</div>';
                        }}
                    }}
                }}
                
                // Start initialization
                initViewer();
            </script>
            
            <div style="margin-top: 15px; padding: 10px; background: #f5f5f5; border-radius: 4px; text-align: center;">
                <div style="font-size: 14px; color: #333; margin-bottom: 5px;">
                    <strong>🧬 {drug_name}</strong> docked to <strong>{target_protein}</strong> protein
                </div>
                <div style="font-size: 12px; color: #666;">
                    🖱️ Rotate: Left-click drag | Zoom: Scroll | Pan: Right-click drag | {len(ligand_data)} poses shown
                </div>
            </div>
        </body>
        </html>
        """
        
        # Render the viewer
        components.html(html_content, height=650, scrolling=False)
        
        return True
        
    except Exception as e:
        logger.error(f"Simple 3D viewer error: {e}")
        st.error(f"3D visualization failed: {e}")
        return False


class Minimal3DViewer:
    """Legacy class for backward compatibility"""
    
    @staticmethod
    def create_3d_view(drug_name: str, target_protein: str):
        """Create 3D view (backward compatible wrapper)"""
        return create_simple_3d_viewer(drug_name, target_protein)
