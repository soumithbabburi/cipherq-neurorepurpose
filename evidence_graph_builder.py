"""
evidence_graph_builder.py - Fixed & Complete Version
Includes get_summary_metrics + rich tooltips for ECharts hover
"""

import pandas as pd
import logging
import os
import json
import sys

logger = logging.getLogger(__name__)

# Global cached data
_CURATED_INTERACTIONS = None
_PATHWAYS = None
_PROTEIN_PATHWAYS = None
_GENES = None


def load_curated_data(force_reload=False):
    """Load all JSON data once"""
    global _CURATED_INTERACTIONS, _PATHWAYS, _PROTEIN_PATHWAYS, _GENES

    if not force_reload and _CURATED_INTERACTIONS is not None:
        return

    base_dir = os.path.dirname(os.path.abspath(__file__))

    try:
        # drug_interactions.json
        path = os.path.join(base_dir, 'drug_interactions.json')
        if os.path.exists(path):
            with open(path, 'r', encoding='utf-8') as f:
                _CURATED_INTERACTIONS = json.load(f)
            logger.info(f"Loaded drug_interactions.json - {len(_CURATED_INTERACTIONS)} entries")
        else:
            logger.warning("drug_interactions.json not found")

        # pathways.json
        path = os.path.join(base_dir, 'pathways.json')
        if os.path.exists(path):
            with open(path, 'r', encoding='utf-8') as f:
                _PATHWAYS = json.load(f)
            logger.info(f"Loaded pathways.json - {len(_PATHWAYS)} pathways")

        # protein_pathways.json
        path = os.path.join(base_dir, 'protein_pathways.json')
        if os.path.exists(path):
            with open(path, 'r', encoding='utf-8') as f:
                _PROTEIN_PATHWAYS = json.load(f)
            logger.info(f"Loaded protein_pathways.json")

        # genes.json
        path = os.path.join(base_dir, 'genes.json')
        if os.path.exists(path):
            with open(path, 'r', encoding='utf-8') as f:
                _GENES = json.load(f)
            logger.info(f"Loaded genes.json - {len(_GENES)} genes")
        else:
            logger.warning("genes.json not found - tooltips will be limited")

    except Exception as e:
        logger.error(f"Error loading JSON data: {e}", exc_info=True)


class EvidenceGraphBuilder:
    def __init__(self):
        load_curated_data()
        logger.info("EvidenceGraphBuilder initialized")

    def build_evidence_graph(self, drug_names: list[str], disease_name: str):
        """
        Main method: build nodes & edges using JSON only
        Returns: (nodes_df, edges_df)
        """
        load_curated_data()

        nodes = []
        edges = []
        node_ids = set()

        # ────────────────────────────────────────────────
        # 1. Add Drug nodes
        # ────────────────────────────────────────────────
        for drug in drug_names:
            drug_id = f"DRUG_{drug.replace(' ', '_').replace('-', '_').upper()}"
            if drug_id not in node_ids:
                tooltip = f"<b>{drug}</b><br>Drug<br>Potential repurposing candidate"
                nodes.append({
                    'id': drug_id,
                    'name': drug,
                    'type': 'drug',
                    'tooltip': tooltip,
                    'symbol': 'circle',
                    'symbolSize': 42
                })
                node_ids.add(drug_id)

        # ────────────────────────────────────────────────
        # 2. Drug → Protein interactions (from drug_interactions.json)
        # ────────────────────────────────────────────────
        interactions = []

        for drug in drug_names:
            key = drug.lower()
            targets = _CURATED_INTERACTIONS.get(key)

            # Flexible matching fallback
            if not targets:
                for k in _CURATED_INTERACTIONS:
                    if key in k.lower() or k.lower() in key:
                        targets = _CURATED_INTERACTIONS[k]
                        logger.info(f"Matched '{drug}' → '{k}'")
                        break

            if targets:
                # Sort by confidence, take top 5-8
                sorted_targets = sorted(
                    targets,
                    key=lambda x: x.get('confidence_score', 0.5),
                    reverse=True
                )[:8]

                for t in sorted_targets:
                    gene = t.get('gene_symbol', '').strip()
                    if not gene:
                        continue

                    protein_name = _GENES.get(gene, {}).get('name', gene)
                    function_snippet = _GENES.get(gene, {}).get('function', '')[:180]
                    if function_snippet:
                        function_snippet += '...' if len(function_snippet) > 170 else ''

                    interactions.append({
                        'drug_name': drug,
                        'gene_symbol': gene,
                        'protein_name': protein_name,
                        'function': function_snippet,
                        'confidence': t.get('confidence_score', 0.7)
                    })

        # ────────────────────────────────────────────────
        # 3. Protein / Gene nodes + edges
        # ────────────────────────────────────────────────
        for inter in interactions:
            gene = inter['gene_symbol']
            prot_id = f"PROTEIN_{gene.upper()}"

            if prot_id not in node_ids:
                tooltip = f"<b>{gene}</b><br>{inter['protein_name']}<br>{inter['function']}<br>Type: Protein target"
                nodes.append({
                    'id': prot_id,
                    'name': gene,
                    'type': 'protein',
                    'tooltip': tooltip,
                    'symbol': 'diamond',
                    'symbolSize': 50
                })
                node_ids.add(prot_id)

            # Drug → Protein edge
            drug_id = f"DRUG_{inter['drug_name'].replace(' ', '_').replace('-', '_').upper()}"
            edges.append({
                'source': drug_id,
                'target': prot_id,
                'label': 'TARGETS',
                'value': inter['confidence']
            })

        # ────────────────────────────────────────────────
        # 4. Pathway nodes & Protein → Pathway edges
        # ────────────────────────────────────────────────
        pathway_nodes_added = set()

        for inter in interactions:
            gene = inter['gene_symbol']
            pathway_ids = _PROTEIN_PATHWAYS.get(gene, [])

            for pid in pathway_ids:
                pinfo = _PATHWAYS.get(pid, {})
                if not pinfo:
                    continue

                pathway_id = f"PATHWAY_{pid}"
                if pathway_id not in node_ids:
                    name = pinfo.get('name', pid)
                    category = pinfo.get('category', 'Unknown')
                    tooltip = f"<b>{name}</b><br>Category: {category}<br>Source: {pinfo.get('source', 'Unknown')}"
                    nodes.append({
                        'id': pathway_id,
                        'name': name[:35] + '...' if len(name) > 35 else name,
                        'full_name': name,
                        'type': 'pathway',
                        'tooltip': tooltip,
                        'symbol': 'rect',
                        'symbolSize': 38
                    })
                    node_ids.add(pathway_id)
                    pathway_nodes_added.add(pathway_id)

                # Protein → Pathway edge
                prot_id = f"PROTEIN_{gene.upper()}"
                edges.append({
                    'source': prot_id,
                    'target': pathway_id,
                    'label': 'PARTICIPATES_IN',
                    'value': 0.85
                })

        # ────────────────────────────────────────────────
        # 5. Disease node + Pathway → Disease connections
        # ────────────────────────────────────────────────
        disease_id = f"DISEASE_{disease_name.replace(' ', '_').replace('-', '_').upper()}"
        if disease_id not in node_ids:
            tooltip = f"<b>{disease_name}</b><br>Target disease"
            nodes.append({
                'id': disease_id,
                'name': disease_name,
                'type': 'disease',
                'tooltip': tooltip,
                'symbol': 'circle',
                'symbolSize': 65
            })
            node_ids.add(disease_id)

        # Connect relevant pathways to disease
        connections = 0
        for pid in pathway_nodes_added:
            pinfo = _PATHWAYS.get(pid.replace("PATHWAY_", ""), {})
            diseases = pinfo.get('diseases', []) or []
            if any(disease_name.lower() in d.lower() for d in diseases):
                edges.append({
                    'source': f"PATHWAY_{pid.replace('PATHWAY_', '')}",
                    'target': disease_id,
                    'label': 'ASSOCIATED_WITH',
                    'value': pinfo.get('relevance_score', 0.7)
                })
                connections += 1

        logger.info(f"Graph built: {len(nodes)} nodes, {len(edges)} edges | {connections} pathway-disease links")

        nodes_df = pd.DataFrame(nodes)
        edges_df = pd.DataFrame(edges)

        return nodes_df, edges_df

    def get_summary_metrics(self, nodes_df: pd.DataFrame, edges_df: pd.DataFrame) -> dict:
        """Restored method - required by your app"""
        if nodes_df.empty:
            return {
                'total_nodes': 0,
                'total_edges': 0,
                'drug_count': 0,
                'protein_count': 0,
                'pathway_count': 0,
                'disease_count': 0,
            }

        type_counts = nodes_df['type'].value_counts().to_dict() if 'type' in nodes_df.columns else {}

        return {
            'total_nodes': len(nodes_df),
            'total_edges': len(edges_df),
            'drug_count': type_counts.get('drug', 0),
            'protein_count': type_counts.get('protein', 0),
            'pathway_count': type_counts.get('pathway', 0),
            'disease_count': type_counts.get('disease', 0),
        }

    def get_pathway_details(self, nodes_df: pd.DataFrame, edges_df: pd.DataFrame):
        """Restored - returns list of pathway info (if your app still uses it)"""
        if nodes_df.empty:
            return []

        pathway_nodes = nodes_df[nodes_df['type'] == 'pathway'].to_dict('records')

        details = []
        for p in pathway_nodes:
            pid = p['id']
            pathway_edges = edges_df[edges_df['target'] == pid]
            proteins = pathway_edges[pathway_edges['label'] == 'PARTICIPATES_IN']['source'].tolist()

            disease_links = edges_df[(edges_df['source'] == pid) & (edges_df['label'].str.contains('ASSOCIATED', na=False))]

            details.append({
                'id': pid,
                'name': p.get('full_name', p.get('name', 'Unknown')),
                'protein_count': len(proteins),
                'proteins': proteins[:5],
                'disease_relevant': len(disease_links) > 0,
                'relevance_score': disease_links['value'].max() if not disease_links.empty else 0.0
            })

        return sorted(details, key=lambda x: x['relevance_score'], reverse=True)
