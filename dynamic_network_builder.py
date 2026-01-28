"""
DYNAMIC NETWORK BUILDER - JSON ONLY VERSION
NO DATABASE - uses drug_interactions.json, genes.json, protein_pathways.json
"""

import logging
import json
from typing import List, Dict, Any, Set

logger = logging.getLogger(__name__)

def get_drug_protein_targets(drug_name: str) -> List[str]:
    """
    Get gene symbols from drug_interactions.json
    """
    try:
        with open('drug_interactions.json', 'r') as f:
            interactions = json.load(f)
        
        # Flexible matching
        drug_lower = drug_name.lower()
        targets = None
        
        # Try exact match
        if drug_lower in interactions:
            targets = interactions[drug_lower]
        else:
            # Try partial match
            for key in interactions.keys():
                if drug_lower in key or key in drug_lower:
                    targets = interactions[key]
                    break
        
        if targets:
            # Return top 5 gene symbols
            sorted_targets = sorted(targets, key=lambda x: x.get('confidence_score', 0), reverse=True)
            gene_symbols = [t['gene_symbol'] for t in sorted_targets[:5]]
            logger.info(f"✅ Targets for {drug_name}: {gene_symbols}")
            return gene_symbols
        else:
            logger.warning(f"⚠️ No targets found for {drug_name}")
            return []
            
    except Exception as e:
        logger.error(f"❌ Query failed for {drug_name}: {e}")
        return []


def create_dynamic_drug_target_network(recommended_drugs: List[Dict], target_disease: str = None) -> Dict:
    """
    Build knowledge graph: DRUG → GENE → PROTEIN → PATHWAY → DISEASE
    Uses JSON files only
    """
    
    if target_disease is None:
        target_disease = 'Alzheimer\'s Disease'
    
    logger.info(f"🧬 Building network for {target_disease} using JSON files")
    
    nodes = []
    links = []
    node_ids = set()
    
    # Process drug names from the recommended list
    drug_names = [drug['name'] if isinstance(drug, dict) else str(drug) for drug in recommended_drugs[:5]]
    
    # 1. DRUG NODES
    for drug_name in drug_names:
        drug_id = f"drug_{drug_name}"
        if drug_id not in node_ids:
            nodes.append({
                "id": drug_id,
                "name": f"{drug_name}\n[DRUG]",
                "category": 0,
                "symbolSize": 65,
                "itemStyle": {"color": "#FF6B6B"},
                "label": {"show": True, "fontSize": 13, "fontWeight": "bold"}
            })
            node_ids.add(drug_id)
    
    # 2. GENE NODES (from JSON)
    gene_set = set()
    drug_gene_links = []
    
    for drug_name in drug_names:
        genes = get_drug_protein_targets(drug_name)
        
        for gene in genes:
            gene_id = f"gene_{gene}"
            if gene_id not in node_ids:
                nodes.append({
                    "id": gene_id,
                    "name": f"{gene}\n[GENE]",
                    "category": 1,
                    "symbolSize": 50,
                    "itemStyle": {"color": "#9B59B6"},
                    "label": {"show": True, "fontSize": 11}
                })
                node_ids.add(gene_id)
                gene_set.add(gene)
            
            drug_gene_links.append({
                "source": f"drug_{drug_name}",
                "target": gene_id,
                "lineStyle": {"width": 2.5, "color": "#9B59B6", "type": "solid"},
                "label": {"show": True, "formatter": "targets"},
                "mechanism": f"{drug_name} targets {gene}"
            })

    # 3. PROTEIN NODES (from genes.json)
    gene_protein_links = []
    protein_set = set()
    
    try:
        with open('genes.json', 'r') as f:
            genes_data = json.load(f)
        
        for gene in gene_set:
            gene_info = genes_data.get(gene.upper(), {})
            protein_name = gene_info.get('name', gene)
            
            protein_id = f"protein_{gene}"
            if protein_id not in node_ids:
                nodes.append({
                    "id": protein_id,
                    "name": f"{protein_name}\n[PROTEIN]",
                    "category": 2,
                    "symbolSize": 50,
                    "itemStyle": {"color": "#4A90E2"},
                    "label": {"show": True, "fontSize": 11}
                })
                node_ids.add(protein_id)
                protein_set.add(gene)
            
            gene_protein_links.append({
                "source": f"gene_{gene}",
                "target": protein_id,
                "lineStyle": {"width": 2, "color": "#4A90E2", "type": "solid"},
                "label": {"show": False},
                "mechanism": f"{gene} encodes {protein_name}"
            })
    except Exception as e:
        logger.warning(f"Could not load genes.json: {e}")
    
    # 4. PATHWAY NODES (based on disease)
    pathways = []
    lower_disease = target_disease.lower()
    if 'alzheimer' in lower_disease or 'dementia' in lower_disease:
        pathways = ['Cholinergic signaling', 'Glutamatergic signaling', 'Neuroinflammation']
    elif 'cardiovascular' in lower_disease or 'hypertension' in lower_disease:
        pathways = ['Renin-Angiotensin system', 'Adrenergic signaling', 'Lipid metabolism']
    elif 'diabetes' in lower_disease:
        pathways = ['Glucose metabolism', 'Insulin signaling', 'PPAR signaling']
    elif 'parkinson' in lower_disease:
        pathways = ['Dopaminergic signaling', 'GABA signaling', 'Neuroinflammation']
    else:
        pathways = ['Therapeutic pathway', 'Drug mechanism']
    
    protein_pathway_links = []
    for pathway in pathways[:3]:
        pathway_id = f"pathway_{pathway}"
        if pathway_id not in node_ids:
            nodes.append({
                "id": pathway_id,
                "name": f"{pathway}\n[PATHWAY]",
                "category": 3,
                "symbolSize": 55,
                "itemStyle": {"color": "#E94B3C"},
                "label": {"show": True, "fontSize": 12, "fontWeight": "bold"}
            })
            node_ids.add(pathway_id)
        
        for protein_name in list(protein_set)[:10]:
            protein_pathway_links.append({
                "source": f"protein_{protein_name}",
                "target": pathway_id,
                "lineStyle": {"width": 2.5, "color": "#E94B3C", "type": "solid"},
                "label": {"show": False},
                "mechanism": f"{protein_name} participates in {pathway}"
            })
    
    # 5. DISEASE NODE
    disease_id = f"disease_{target_disease}"
    nodes.append({
        "id": disease_id,
        "name": f"{target_disease}\n[DISEASE]",
        "category": 4,
        "symbolSize": 75,
        "itemStyle": {"color": "#FF0000"},
        "label": {"show": True, "fontSize": 15, "fontWeight": "bold"}
    })
    
    # 6. PATHWAY → DISEASE LINKS
    pathway_disease_links = []
    for pathway in pathways:
        pathway_disease_links.append({
            "source": f"pathway_{pathway}",
            "target": disease_id,
            "lineStyle": {"width": 3.5, "color": "#2ECC71", "type": "solid"},
            "label": {"show": True, "formatter": "dysregulated in", "fontSize": 10},
            "mechanism": f"Pathway dysfunction in {target_disease}"
        })
    
    # 7. ASSEMBLE ALL LINKS
    links.extend(drug_gene_links)
    links.extend(gene_protein_links)
    links.extend(protein_pathway_links)
    links.extend(pathway_disease_links)
    
    logger.info(f"✅ Built network: {len(nodes)} nodes, {len(links)} links")
    
    return {
        "graph_data": {"nodes": nodes, "links": links},
        "title": f"Drug-Gene-Protein Network: {target_disease}",
        "disease": target_disease,
        "summary": f"{len(drug_names)} drugs → {len(gene_set)} genes → {len(protein_set)} proteins"
    }

__all__ = [
    'create_dynamic_drug_target_network',
    'get_drug_protein_targets'
]
