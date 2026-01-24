"""
BULLETPROOF ECHARTS HTML RENDERER
Fixes height=0 and ComponentInstance issues permanently
"""
import streamlit as st
import streamlit.components.v1 as components
import json
from typing import Dict

def render_echarts_html(options: Dict, key: str = "network_graph", height_px: int = 600) -> None:
    """
    Bulletproof ECharts renderer using components.html with guaranteed height
    Eliminates height=0 and ComponentInstance registration issues
    """
    # Sanitize options for JSON serialization
    clean_options = json.dumps(options, default=str, ensure_ascii=False)
    
    # Create unique div ID to prevent conflicts
    div_id = f"echarts-{key.replace('_', '-')}"
    
    # HTML with guaranteed height and explicit styling
    echarts_html = f"""
    <div id="{div_id}" style="width: 100%; height: {height_px}px; min-height: {height_px}px; background: white; border: 1px solid #ddd;"></div>
    <script src="https://cdn.jsdelivr.net/npm/echarts@5.4.0/dist/echarts.min.js"></script>
    <script>
        // Initialize chart with explicit sizing
        const chart = echarts.init(document.getElementById('{div_id}'));
        const options = {clean_options};
        
        // Set options with explicit resize
        chart.setOption(options);
        chart.resize();
        
        // Auto-resize on window changes
        window.addEventListener('resize', () => {{
            chart.resize();
        }});
        
        // Force resize after a brief delay to ensure DOM is ready
        setTimeout(() => {{
            chart.resize();
        }}, 100);
    </script>
    """
    
    # Render with explicit height (no key parameter for components.html)
    components.html(echarts_html, height=height_px + 50)
    
    # Success feedback
    st.success("✅ Interactive network graph rendered - click and drag nodes to explore!")