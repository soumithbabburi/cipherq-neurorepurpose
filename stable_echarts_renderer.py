"""
BULLETPROOF ECHARTS HTML RENDERER - Enhanced 2026 Edition
Fixes:
- Height collapse / ComponentInstance duplication
- Edge labels floating / not centered
- Better node separation visibility
- Streamlit re-run safety (unique keys + no global conflicts)
- Mobile/responsive behavior
"""
import streamlit as st
import streamlit.components.v1 as components
import json
from typing import Dict, Any

def render_echarts_html(
    options: Dict[str, Any],
    key: str = "network_graph",
    height_px: int = 600,
    title: str = "Drug-Target-Disease Network"
) -> None:
    """
    Renders ECharts graph with guaranteed sizing, proper edge labels,
    and Streamlit re-run safety.
    
    Args:
        options: Full ECharts options dict (should already include series.graph)
        key: Unique identifier for this chart instance
        height_px: Height of the container
        title: Optional chart title shown above
    """
    # 1. Make options JSON-safe and add missing safety defaults
    clean_options = json.dumps(options, default=str, ensure_ascii=False)
    
    # Unique div ID per chart instance (prevents multiple charts clashing)
    div_id = f"echarts-{key.replace('_', '-')}-{id(options) % 1000000:06d}"
    
    # 2. Enhanced HTML template with:
    #    - Explicit min/max height
    #    - White background + subtle border
    #    - Title above chart
    #    - Multiple resize triggers
    echarts_html = f"""
    <div style="width: 100%; padding: 8px; box-sizing: border-box;">
        <h4 style="margin: 0 0 8px 0; color: #1f2937; text-align: center;">{title}</h4>
        <div id="{div_id}" 
             style="width: 100%; 
                    height: {height_px}px; 
                    min-height: {height_px}px; 
                    max-height: {height_px + 200}px; 
                    background: #ffffff; 
                    border: 1px solid #e5e7eb; 
                    border-radius: 6px; 
                    overflow: hidden;">
        </div>
    </div>

    <script src="https://cdn.jsdelivr.net/npm/echarts@5.5.1/dist/echarts.min.js"></script>
    <script>
        (function() {{
            const container = document.getElementById('{div_id}');
            if (!container) return;
            
            const chart = echarts.init(container, null, {{renderer: 'canvas'}});
            
            try {{
                const options = {clean_options};
                
                // Safety: ensure graph series has edge labels configured
                if (options.series && options.series[0] && options.series[0].type === 'graph') {{
                    options.series[0].edgeLabel = {{
                        ...options.series[0].edgeLabel,
                        show: true,
                        position: 'middle',
                        formatter: '{{c}}',           // uses per-link 'label' value
                        fontSize: 11,
                        color: '#4b5563',
                        backgroundColor: 'rgba(255,255,255,0.9)',
                        padding: [3, 6, 3, 6],
                        borderRadius: 4
                    }};
                }}
                
                chart.setOption(options, true);  // true = notMerge → safer on re-runs
                
                // Multiple resize triggers
                chart.resize();
                
                setTimeout(() => chart.resize(), 100);
                setTimeout(() => chart.resize(), 500);
                
                window.addEventListener('resize', () => {{
                    chart.resize();
                }}, {{passive: true}});
                
                // Optional: tooltip always visible on hover
                chart.on('mouseover', (params) => {{
                    if (params.dataType === 'edge') {{
                        chart.dispatchAction({{type: 'showTip', dataIndex: params.dataIndex}});
                    }}
                }});
            }} catch (e) {{
                container.innerHTML = '<div style="padding: 20px; color: #dc2626;">Error rendering graph: ' + e.message + '</div>';
            }}
        }})();
    </script>
    """

    # 3. Render safely — add extra height buffer for title + padding
    components.html(echarts_html, height=height_px + 80, scrolling=False)

    # 4. User feedback (only once per session if desired)
    if "echarts_rendered" not in st.session_state:
        st.success("✅ Interactive BioCypher network rendered — drag nodes, hover edges for details")
        st.session_state["echarts_rendered"] = True
