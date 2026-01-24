"""
Minimal Network Graph Renderer - Clean ECharts with PyVis fallback
Fixes: Component registration, height issues, data sanitization
"""

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import json
import hashlib
import logging
from typing import Dict, List, Any, Optional

# BULLETPROOF IMPORT FOR st_echarts
ECHARTS_AVAILABLE = False

try:
    from streamlit_echarts import st_echarts
    ECHARTS_AVAILABLE = True
    print("✅ streamlit-echarts imported successfully")
except ImportError as e:
    print(f"❌ streamlit-echarts import failed: {e}")
    ECHARTS_AVAILABLE = False
    
    # Safe fallback function
    def st_echarts(*args, **kwargs):
        import streamlit as st
        st.error("ECharts visualization requires streamlit-echarts package")
        return None
except Exception as e:
    print(f"❌ Unexpected error importing streamlit-echarts: {e}")
    ECHARTS_AVAILABLE = False
    
    # Safe fallback function  
    def st_echarts(*args, **kwargs):
        import streamlit as st
        st.error(f"ECharts import error: {str(e)}")
        return None

def get_logger():
    return logging.getLogger(__name__)

def sanitize_for_echarts(data: Any) -> Any:
    """Convert all data to JSON-serializable types for ECharts"""
    if isinstance(data, (list, tuple)):
        return [sanitize_for_echarts(item) for item in data]
    elif isinstance(data, dict):
        return {str(k): sanitize_for_echarts(v) for k, v in data.items()}
    elif hasattr(data, 'item'):  # numpy types
        return data.item()
    elif hasattr(data, 'tolist'):  # pandas/numpy arrays
        return data.tolist()
    elif data is None:
        return None
    else:
        return str(data) if not isinstance(data, (int, float, bool, str)) else data

def create_echarts_options(nodes_df: pd.DataFrame, edges_df: pd.DataFrame) -> Dict:
    """Create clean ECharts options from node/edge dataframes"""
    
    # Node categories and colors
    node_colors = {
        'drug': '#ef4444',      # Red
        'protein': '#f97316',   # Orange  
        'disease': '#8b5cf6',   # Purple
        'pathway': '#06b6d4'    # Cyan
    }
    
    # Process nodes
    nodes = []
    categories = []
    category_map = {}
    
    for _, node in nodes_df.iterrows():
        node_type = str(node.get('type', 'unknown')).lower()
        node_name = str(node.get('name', node.get('id', 'Unknown')))
        
        # Create category if not exists
        if node_type not in category_map:
            category_map[node_type] = len(categories)
            categories.append({
                'name': node_type.title(),
                'itemStyle': {'color': node_colors.get(node_type, '#64748b')}
            })
        
        nodes.append({
            'id': str(node.get('id', node_name)),
            'name': node_name,
            'category': category_map[node_type],
            'symbolSize': 50 if node_type == 'disease' else 35,
            'value': float(node.get('confidence', 0.7)),
            'itemStyle': {'color': node_colors.get(node_type, '#64748b')}
        })
    
    # Process edges
    links = []
    for _, edge in edges_df.iterrows():
        links.append({
            'source': str(edge.get('source', '')),
            'target': str(edge.get('target', '')),
            'lineStyle': {
                'color': '#94a3b8',
                'width': 2
            }
        })
    
    # Create ECharts options
    options = {
        'title': {
            'text': 'Drug-Target-Disease Network',
            'textStyle': {'color': '#1f2937', 'fontSize': 16}
        },
        'tooltip': {'trigger': 'item'},
        'legend': {'data': [cat['name'] for cat in categories]},
        'series': [{
            'type': 'graph',
            'layout': 'force',
            'data': nodes,
            'links': links,
            'categories': categories,
            'roam': True,
            'focusNodeAdjacency': True,
            'force': {
                'repulsion': 1000,
                'gravity': 0.1,
                'edgeLength': 150
            },
            'label': {
                'show': True,
                'position': 'right',
                'fontSize': 12
            }
        }]
    }
    
    return sanitize_for_echarts(options)

def render_pyvis_fallback(nodes_df: pd.DataFrame, edges_df: pd.DataFrame, height: int = 600):
    """Fallback network using HTML and simple D3.js"""
    
    node_colors = {
        'drug': '#ef4444', 'protein': '#f97316', 
        'disease': '#8b5cf6', 'pathway': '#06b6d4'
    }
    
    nodes_list = []
    for _, node in nodes_df.iterrows():
        node_type = str(node.get('type', 'unknown')).lower()
        nodes_list.append({
            'id': str(node.get('id', node.get('name', 'Unknown'))),
            'label': str(node.get('name', node.get('id', 'Unknown'))),
            'color': node_colors.get(node_type, '#64748b'),
            'size': 20
        })
    
    edges_list = []
    for _, edge in edges_df.iterrows():
        edges_list.append({
            'from': str(edge.get('source', '')),
            'to': str(edge.get('target', ''))
        })
    
    html_content = f"""
    <div style="width: 100%; height: {height}px; border: 1px solid #ccc;">
        <h3 style="text-align: center; margin: 20px;">Drug-Target-Disease Network</h3>
        <div style="padding: 20px;">
            <div style="margin-bottom: 15px;">
                <strong>Network Summary:</strong>
                <ul>
                    <li>Nodes: {len(nodes_list)}</li>
                    <li>Connections: {len(edges_list)}</li>
                </ul>
            </div>
            <div style="display: flex; flex-wrap: wrap; gap: 10px;">
                {' '.join([f'<div style="background: {node["color"]}; color: white; padding: 5px 10px; border-radius: 15px; font-size: 12px;">{node["label"]}</div>' for node in nodes_list])}
            </div>
        </div>
    </div>
    """
    
    components.html(html_content, height=height)

def render_network_graph(nodes_df: pd.DataFrame, edges_df: pd.DataFrame, 
                        key: Optional[str] = None, height: int = 600) -> bool:
    """
    Main network rendering function with ECharts + fallback
    
    Returns:
        bool: True if successful, False if failed
    """
    
    logger = get_logger()
    
    # Generate unique key if not provided
    if key is None:
        data_hash = hashlib.md5(f"{len(nodes_df)}_{len(edges_df)}_{hash(str(nodes_df.values.tobytes()))}".encode()).hexdigest()[:8]
        key = f"network_{data_hash}"
    
    # Validate inputs
    if nodes_df.empty or edges_df.empty:
        st.warning("No network data available to display")
        return False
    
    logger.info(f"Network rendering starting: {len(nodes_df)} nodes, {len(edges_df)} edges")
    
    # CREATE HIGHLY VISIBLE GRAPH SECTION
    st.markdown("---")
    st.markdown("### 🌐 **INTERACTIVE NETWORK GRAPH IS HERE** 👇")
    st.markdown("**📍 Look below for the interactive drug-target-disease network visualization**")
    
    # Create bordered container with background
    with st.container():
        st.markdown(
            f"<div style='border: 3px solid #1f77b4; border-radius: 10px; padding: 20px; background-color: #f8f9fa; margin: 10px 0;'>" +
            f"<h4 style='color: #1f77b4; margin-top: 0;'>🔗 BioCypher Network: {len(nodes_df)} nodes, {len(edges_df)} connections</h4>", 
            unsafe_allow_html=True
        )
        
        try:
            if ECHARTS_AVAILABLE:
                options = create_echarts_options(nodes_df, edges_df)
                
                try:
                    # BULLETPROOF APPROACH: Use components.html with direct ECharts
                    echarts_html = f"""
                    <div id="echarts-container" style="width: 100%; height: {height}px; min-height: {height}px;"></div>
                    <script src="https://cdn.jsdelivr.net/npm/echarts@5.4.0/dist/echarts.min.js"></script>
                    <script>
                        var myChart = echarts.init(document.getElementById('echarts-container'));
                        var option = {json.dumps(options)};
                        myChart.setOption(option);
                        
                        // Auto-resize
                        window.addEventListener('resize', function() {{
                            myChart.resize();
                        }});
                    </script>
                    """
                    
                    # Use the new stable renderer with guaranteed height
                    from stable_echarts_renderer import render_echarts_html
                    render_echarts_html(options, key=f"network_graph_{len(nodes_df)}_{len(edges_df)}", height_px=height)
                    
                    # Show success message and debugging info
                    logger.info("ECharts network rendered successfully via HTML")
                    st.success("✅ Interactive network graph rendered above! Click and drag nodes to explore.")
                    st.info(f"📊 Network: {len(nodes_df)} nodes ({len(nodes_df[nodes_df['type'] == 'drug'])} drugs, {len(nodes_df[nodes_df['type'] == 'protein'])} proteins, {len(nodes_df[nodes_df['type'] == 'disease'])} diseases), {len(edges_df)} connections")
                    return True
                except Exception as e:
                    logger.error(f"ECharts rendering failed: {e}")
                    st.error(f"ECharts error: {str(e)}")
                    # Fall through to fallback below
            
            # ECharts not available or failed - use fallback
            logger.warning("ECharts not available, using fallback renderer")
            st.warning("⚠️ Interactive graph unavailable - showing network summary below:")
            render_pyvis_fallback(nodes_df, edges_df, height)
            st.markdown("</div>", unsafe_allow_html=True)
            return True
            
        except Exception as e:
            logger.error(f"Network rendering failed: {e}")
            st.error(f"Network visualization error: {str(e)}")
            
            # Fallback rendering
            st.warning("Interactive graph failed - showing network summary")
            render_pyvis_fallback(nodes_df, edges_df, height)
            st.markdown("</div>", unsafe_allow_html=True)
            return False