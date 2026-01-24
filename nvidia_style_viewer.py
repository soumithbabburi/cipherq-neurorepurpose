"""
Professional NVIDIA DiffDock-Style 3D Viewer
Beautiful rendering with pose ranking sidebar
"""
import streamlit as st
import streamlit.components.v1 as components
import logging
import os
import glob

logger = logging.getLogger(__name__)

def create_nvidia_style_viewer(drug_name: str, target_protein: str, poses: list, height: int = 700):
    """
    Create professional DiffDock-style viewer with:
    - Dark background
    - Multiple poses shown
    - Ranking sidebar
    - Professional protein rendering
    """
    
    try:
        # Get protein PDB
        diffdock_dir = f"./diffdock_output/{drug_name}"
        protein_pdb_path = f"{diffdock_dir}/protein.pdb"
        
        if not os.path.exists(protein_pdb_path):
            # Try pdb_cache
            pdb_files = glob.glob(f"./pdb_cache/*.pdb")
            if pdb_files:
                protein_pdb_path = pdb_files[0]
            else:
                st.error(f"No protein structure found for {target_protein}")
                return False
        
        # Read protein PDB
        with open(protein_pdb_path, 'r') as f:
            protein_pdb = f.read()
        
        # Get ligand SDF files (up to 5 poses for visualization)
        ligand_files = glob.glob(f"{diffdock_dir}/pose_*.sdf")[:5]
        
        if not ligand_files:
            st.warning(f"No ligand poses found for {drug_name}")
            return False
        
        # Read all ligand poses
        ligand_sdfs = []
        for sdf_file in ligand_files:
            with open(sdf_file, 'r') as f:
                ligand_sdfs.append(f.read())
        
        # Escape for JavaScript
        protein_pdb_js = protein_pdb.replace('`', '\\`').replace('${', '\\${')
        ligand_sdfs_js = [sdf.replace('`', '\\`').replace('${', '\\${') for sdf in ligand_sdfs]
        
        # Build pose ranking list
        pose_list_html = ""
        for i, pose in enumerate(poses[:20]):
            affinity = pose.get('binding_affinity', 0)
            score = pose.get('confidence', 0)
            pose_list_html += f"""
            <div class="pose-item" onclick="showPose({i})">
                <span class="rank">Rank: {i+1}</span>
                <span class="score">Score: {score:.3f}</span>
            </div>
            """
        
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
            <style>
                * {{ margin: 0; padding: 0; box-sizing: border-box; }}
                body {{ 
                    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
                    background: #000;
                    color: #fff;
                }}
                .container {{
                    display: flex;
                    height: {height}px;
                    background: #000;
                }}
                #viewer {{
                    flex: 1;
                    position: relative;
                }}
                .sidebar {{
                    width: 280px;
                    background: #1a1a1a;
                    border-left: 1px solid #333;
                    overflow-y: auto;
                    padding: 20px;
                }}
                .sidebar h3 {{
                    font-size: 18px;
                    margin-bottom: 15px;
                    color: #9ef01a;
                    display: flex;
                    align-items: center;
                    gap: 8px;
                }}
                .pose-item {{
                    background: #2a2a2a;
                    padding: 12px;
                    margin-bottom: 8px;
                    border-radius: 6px;
                    cursor: pointer;
                    transition: all 0.2s;
                    border: 1px solid #333;
                }}
                .pose-item:hover {{
                    background: #3a3a3a;
                    border-color: #9ef01a;
                }}
                .pose-item.active {{
                    background: #3a3a3a;
                    border-color: #9ef01a;
                }}
                .rank {{
                    display: block;
                    font-size: 14px;
                    color: #aaa;
                }}
                .score {{
                    display: block;
                    font-size: 16px;
                    font-weight: 600;
                    color: #fff;
                    margin-top: 4px;
                }}
                .info-overlay {{
                    position: absolute;
                    top: 20px;
                    left: 20px;
                    background: rgba(26, 26, 26, 0.95);
                    padding: 15px 20px;
                    border-radius: 8px;
                    border: 1px solid #333;
                    z-index: 1000;
                }}
                .info-title {{
                    font-size: 16px;
                    font-weight: 600;
                    color: #9ef01a;
                    margin-bottom: 8px;
                }}
                .info-item {{
                    font-size: 13px;
                    color: #aaa;
                    margin: 4px 0;
                }}
                .controls {{
                    position: absolute;
                    bottom: 20px;
                    left: 20px;
                    background: rgba(26, 26, 26, 0.9);
                    padding: 10px 15px;
                    border-radius: 6px;
                    font-size: 12px;
                    color: #aaa;
                }}
            </style>
        </head>
        <body>
            <div class="container">
                <div id="viewer">
                    <div class="info-overlay">
                        <div class="info-title">🧬 {drug_name} → {target_protein}</div>
                        <div class="info-item" id="current-pose">Pose: #1</div>
                        <div class="info-item" id="current-affinity">Affinity: {poses[0].get('binding_affinity', 0):.1f} kcal/mol</div>
                    </div>
                    <div class="controls">
                        🖱️ Rotate: Left-drag | Zoom: Scroll | Pan: Right-drag
                    </div>
                </div>
                <div class="sidebar">
                    <h3>
                        <svg width="20" height="20" viewBox="0 0 20 20" fill="#9ef01a">
                            <circle cx="10" cy="10" r="8"/>
                        </svg>
                        View All Poses
                    </h3>
                    {pose_list_html}
                </div>
            </div>
            
            <script>
                let viewer = $3Dmol.createViewer(document.getElementById('viewer'), {{
                    backgroundColor: '#000000'
                }});
                
                // Load protein (only once)
                let proteinData = `{protein_pdb_js}`;
                viewer.addModel(proteinData, "pdb");
                
                // Style protein as ribbons (NVIDIA style)
                viewer.setStyle({{model: 0}}, {{
                    cartoon: {{
                        color: 'lightgray',
                        opacity: 0.8,
                        thickness: 0.5
                    }}
                }});
                
                // Add all ligand poses
                let ligands = [
                    {', '.join([f'`{sdf}`' for sdf in ligand_sdfs_js])}
                ];
                
                ligands.forEach((ligandData, idx) => {{
                    viewer.addModel(ligandData, "sdf");
                    // Initially hide all except first
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
                
                // Center on ligands
                viewer.zoomTo({{model: 1}});
                viewer.zoom(0.8, 500);
                viewer.render();
                
                // Function to switch poses
                window.showPose = function(poseIndex) {{
                    // Hide all ligands
                    for (let i = 0; i < ligands.length; i++) {{
                        viewer.setStyle({{model: i + 1}}, {{
                            stick: {{ hidden: true }},
                            sphere: {{ hidden: true }}
                        }});
                    }}
                    
                    // Show selected ligand
                    viewer.setStyle({{model: poseIndex + 1}}, {{
                        stick: {{
                            colorscheme: 'Jmol',
                            radius: 0.25,
                            hidden: false
                        }},
                        sphere: {{
                            scale: 0.35,
                            colorscheme: 'Jmol',
                            hidden: false
                        }}
                    }});
                    
                    viewer.render();
                    
                    // Update info overlay
                    let poses = {[f'{{"affinity": {p.get("binding_affinity", 0)}}}' for p in poses[:20]]};
                    document.getElementById('current-pose').textContent = 'Pose: #' + (poseIndex + 1);
                    document.getElementById('current-affinity').textContent = 'Affinity: ' + poses[poseIndex].affinity.toFixed(1) + ' kcal/mol';
                    
                    // Highlight active pose in sidebar
                    document.querySelectorAll('.pose-item').forEach((item, idx) => {{
                        if (idx === poseIndex) {{
                            item.classList.add('active');
                        }} else {{
                            item.classList.remove('active');
                        }}
                    }});
                }};
                
                // Highlight first pose
                document.querySelector('.pose-item').classList.add('active');
            </script>
        </body>
        </html>
        """
        
        components.html(html, height=height)
        
        st.caption(f"💡 **Click poses on the right to switch between different binding configurations**")
        
        logger.info(f"✅ Professional NVIDIA-style viewer rendered: {drug_name} → {target_protein}")
        return True
        
    except Exception as e:
        logger.error(f"Professional viewer failed: {e}")
        st.error(f"3D visualization error: {e}")
        return False
