import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import json
from stable_echarts_renderer import render_echarts_html

def create_echarts_options(nodes_df: pd.DataFrame, edges_df: pd.DataFrame):
    node_colors = {
        'drug': '#3b82f6',
        'protein': '#10b981',   # Green for proteins/genes
        'pathway': '#f59e0b',
        'disease': '#ef4444'
    }

    nodes = []
    for _, row in nodes_df.iterrows():
        ntype = str(row.get('type', 'unknown')).lower()
        nodes.append({
            'id': str(row.get('id')),
            'name': str(row.get('name')),
            'category': ntype.title(),
            'symbol': row.get('symbol', 'circle'),
            'symbolSize': int(row.get('symbolSize', 35)),
            'itemStyle': {'color': node_colors.get(ntype, '#64748b')},
            'tooltip': row.get('tooltip', str(row.get('name')))  # ← Rich tooltip
        })

    links = [{'source': str(e.get('source')), 'target': str(e.get('target'))} for _, e in edges_df.iterrows()]

    options = {
        'title': {'text': 'Drug-Target-Disease Network', 'textStyle': {'fontSize': 18}},
        'tooltip': {
            'trigger': 'item',
            'formatter': '''
                function(params) {
                    return params.data.tooltip || params.name;
                }
            '''
        },
        'legend': {'data': list(node_colors.keys())},
        'series': [{
            'type': 'graph',
            'layout': 'force',
            'data': nodes,
            'links': links,
            'categories': [{'name': k.title()} for k in node_colors],
            'roam': True,
            'focusNodeAdjacency': True,
            'force': {
                'repulsion': 2200,
                'gravity': 0.05,
                'edgeLength': [120, 280]
            },
            'label': {
                'show': True,
                'position': 'right',
                'fontSize': 12,
                'color': '#1f2937',
                'textBorderColor': '#ffffff',
                'textBorderWidth': 3
            },
            'emphasis': {'focus': 'adjacency', 'lineStyle': {'width': 4}}
        }]
    }
    return options

def render_network_graph(nodes_df, edges_df, height=700):
    if nodes_df.empty:
        st.warning("No data")
        return
    st.markdown("### Interactive Drug-Target-Disease Network")
    options = create_echarts_options(nodes_df, edges_df)
    render_echarts_html(options, key="fixed_network", height_px=height)
    st.caption("💡 Drag nodes, zoom with mouse wheel, hover for descriptions")
