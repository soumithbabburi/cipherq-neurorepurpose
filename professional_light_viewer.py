"""
Professional 3D Viewer - NVIDIA Quality with Clean Light Theme
Full protein ribbons, multiple poses, BUT with your clean white interface
"""
import streamlit as st
import streamlit.components.v1 as components
import logging
import os
import glob

logger = logging.getLogger(__name__)

def create_professional_light_viewer(drug_name: str, target_protein: str, poses: list, height: int = 650):
    """
    Professional 3D viewer:
    - NVIDIA-quality protein rendering (full ribbons, helices, sheets)
    - Light/white background (matches your interface)
    - Multiple poses switchable
    - Professional appearance
    """
    
    try:
        # Get protein PDB
        diffdock_dir = f"./diffdock_output/{drug_name}"
        protein_pdb_path = f"{diffdock_dir}/protein.pdb"
        
        if not os.path.exists(protein_pdb_path):
            pdb_files = glob.glob(f"./pdb_cache/*.pdb")
            if pdb_files:
                protein_pdb_path = pdb_files[0]
            else:
                st.warning(f"No protein structure available for {target_protein}")
                return False
        
        with open(protein_pdb_path, 'r') as f:
            protein_pdb = f.read()
        
        # Get ligand poses (up to 10)
        ligand_files = glob.glob(f"{diffdock_dir}/pose_*.sdf")[:10]
        
        if not ligand_files:
            st.warning(f"No ligand poses found")
            return False
        
        ligand_sdfs = []
        for sdf_file in ligand_files:
            with open(sdf_file, 'r') as f:
                ligand_sdfs.append(f.read())
        
        # Escape for JavaScript
        protein_pdb_js = protein_pdb.replace('`', '\\`').replace('${', '\\${')
        ligand_sdfs_js = [sdf.replace('`', '\\`').replace('${', '\\${') for sdf in ligand_sdfs]
        
        # Build pose buttons (horizontal)
        pose_buttons = ""
        for i, pose in enumerate(poses[:10]):
            affinity = pose.get('binding_affinity', 0)
            pose_buttons += f'<button class="pose-btn" onclick="showPose({i})">Pose {i+1}<br>{affinity:.1f} kcal/mol</button>'
        
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
            <style>
                body {{ margin: 0; padding: 0; background: white; }}
                #viewer-container {{ width: 100%; height: {height}px; position: relative; background: white; }}
                #viewer {{ width: 100%; height: 100%; }}
                .info-box {{
                    position: absolute;
                    top: 15px;
                    left: 15px;
                    background: rgba(255, 255, 255, 0.95);
                    padding: 12px 18px;
                    border-radius: 8px;
                    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
                    border: 1px solid #e2e8f0;
                    font-family: -apple-system, sans-serif;
                    z-index: 1000;
                }}
                .info-title {{ font-size: 15px; font-weight: 600; color: #1e293b; margin-bottom: 6px; }}
                .info-detail {{ font-size: 13px; color: #64748b; margin: 3px 0; }}
                .pose-controls {{
                    position: absolute;
                    bottom: 15px;
                    left: 50%;
                    transform: translateX(-50%);
                    background: rgba(255, 255, 255, 0.95);
                    padding: 12px;
                    border-radius: 8px;
                    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
                    border: 1px solid #e2e8f0;
                    display: flex;
                    gap: 8px;
                    z-index: 1000;
                }}
                .pose-btn {{
                    background: #f8fafc;
                    border: 1px solid #cbd5e1;
                    padding: 8px 12px;
                    border-radius: 6px;
                    cursor: pointer;
                    font-size: 12px;
                    line-height: 1.3;
                    text-align: center;
                    transition: all 0.2s;
                }}
                .pose-btn:hover {{ background: #e0e7ff; border-color: #6366f1; }}
                .pose-btn.active {{ background: #6366f1; color: white; border-color: #6366f1; font-weight: 600; }}
            </style>
        </head>
        <body>
            <div id="viewer-container">
                <div class="info-box">
                    <div class="info-title">🧬 {drug_name} → {target_protein}</div>
                    <div class="info-detail" id="current-pose">Pose: #1</div>
                    <div class="info-detail" id="current-affinity">Affinity: {poses[0].get('binding_affinity', 0):.1f} kcal/mol</div>
                </div>
                <div id="viewer"></div>
                <div class="pose-controls">
                    {pose_buttons}
                </div>
            </div>
            
            <script>
                let viewer = $3Dmol.createViewer(document.getElementById('viewer'), {{
                    backgroundColor: 'white'  // LIGHT background!
                }});
                
                // Load protein with FULL structure rendering
                let proteinData = `{protein_pdb_js}`;
                viewer.addModel(proteinData, "pdb");
                
                // NVIDIA-QUALITY protein rendering (full ribbons, helices, sheets)
                viewer.setStyle({{model: 0}}, {{
                    cartoon: {{
                        color: 'spectrum',  // Rainbow colors show structure
                        opacity: 0.8,
                        thickness: 1.0,  // Thick ribbons
                        arrows: true,  // Show beta sheets as arrows
                        tubes: true  // Show helices as tubes
                    }}
                }});
                
                // Add surface for context
                viewer.addSurface($3Dmol.SurfaceType.VDW, {{
                    opacity: 0.1,
                    color: 'lightgray'
                }}, {{model: 0}});
                
                // Load all ligand poses
                let ligands = [{', '.join([f'`{sdf}`' for sdf in ligand_sdfs_js])}];
                
                ligands.forEach((ligandData, idx) => {{
                    viewer.addModel(ligandData, "sdf");
                    viewer.setStyle({{model: idx + 1}}, {{
                        stick: {{
                            colorscheme: 'Jmol',
                            radius: 0.25,
                            hidden: idx !== 0
                        }},
                        sphere: {{
                            scale: 0.35,
                            colorscheme: 'Jmol',
                            hidden: idx !== 0
                        }}
                    }});
                }});
                
                // Center on binding site
                viewer.zoomTo({{model: 1}});
                viewer.zoom(0.75, 500);
                viewer.render();
                
                // Pose switching function
                window.showPose = function(poseIndex) {{
                    // Hide all ligands
                    for (let i = 0; i < ligands.length; i++) {{
                        viewer.setStyle({{model: i + 1}}, {{ stick: {{ hidden: true }}, sphere: {{ hidden: true }} }});
                    }}
                    
                    // Show selected
                    viewer.setStyle({{model: poseIndex + 1}}, {{
                        stick: {{ colorscheme: 'Jmol', radius: 0.25, hidden: false }},
                        sphere: {{ scale: 0.35, colorscheme: 'Jmol', hidden: false }}
                    }});
                    
                    viewer.render();
                    
                    // Update info
                    let poses = {[{"affinity": pose.get("binding_affinity", 0)} for pose in poses[:10]]};
                    document.getElementById('current-pose').textContent = 'Pose: #' + (poseIndex + 1);
                    document.getElementById('current-affinity').textContent = 'Affinity: ' + poses[poseIndex].affinity.toFixed(1) + ' kcal/mol';
                    
                    // Highlight active button
                    document.querySelectorAll('.pose-btn').forEach((btn, idx) => {{
                        if (idx === poseIndex) {{
                            btn.classList.add('active');
                        }} else {{
                            btn.classList.remove('active');
                        }}
                    }});
                }};
                
                // Activate first pose
                document.querySelector('.pose-btn').classList.add('active');
            </script>
        </body>
        </html>
        """
        
        components.html(html, height=height + 50)
        
        logger.info(f"✅ Professional viewer: {drug_name} → {target_protein}")
        return True
        
    except Exception as e:
        logger.error(f"Professional viewer error: {e}")
        return False
