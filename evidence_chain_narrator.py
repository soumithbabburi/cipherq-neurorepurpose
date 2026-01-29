"""
EVIDENCE CHAIN NARRATOR
Generates narrative description of Drug → Protein → Pathway → Disease connections
"""

import pandas as pd
import logging

logger = logging.getLogger(__name__)

def generate_evidence_chain_description(nodes_df, edges_df, drug_name: str, disease_name: str) -> str:
    """
    Generate narrative description of complete evidence chain
    Uses genes.json for full protein names - NO HARDCODING
    """
    
    try:
        import json
        
        # Load genes.json for full protein names
        try:
            with open('genes.json', 'r') as f:
                genes_data = json.load(f)
        except:
            genes_data = {}
        
        # Find drug node
        drug_id = f"DRUG_{drug_name.replace(' ', '_').replace('-', '_').upper()}"
        
        # Get targets (drug → protein edges)
        drug_edges = edges_df[edges_df['source'] == drug_id]
        
        if len(drug_edges) == 0:
            return f"No evidence chain found for {drug_name}."
        
        # Get primary target (first/highest confidence)
        primary_target_id = drug_edges.iloc[0]['target']
        target_node = nodes_df[nodes_df['id'] == primary_target_id]
        
        if len(target_node) == 0:
            return f"{drug_name} targets unknown proteins."
        
        target_gene_symbol = target_node.iloc[0]['name']
        
        # Get FULL protein name from genes.json
        gene_info = genes_data.get(target_gene_symbol.upper(), {})
        full_protein_name = gene_info.get('name', target_gene_symbol)
        
        # Get pathways (protein → pathway edges)
        protein_edges = edges_df[edges_df['source'] == primary_target_id]
        
        pathways = []
        for _, edge in protein_edges.iterrows():
            pathway_id = edge['target']
            pathway_node = nodes_df[nodes_df['id'] == pathway_id]
            if len(pathway_node) > 0 and pathway_node.iloc[0].get('type') == 'pathway':
                pathways.append(pathway_node.iloc[0]['name'])
        
        # Get disease connections
        disease_id = f"DISEASE_{disease_name.replace(' ', '_').replace('-', '_').upper()}"
        pathway_edges_to_disease = edges_df[edges_df['target'] == disease_id]
        connected_pathway_count = len(pathway_edges_to_disease)
        
        # Build narrative
        description = f"**{drug_name} → {disease_name} Evidence Chain:**\n\n"
        
        # Drug → Protein (using full name from genes.json!)
        description += f"🔹 **{drug_name}** TARGETS **{target_gene_symbol}** ({full_protein_name})\n\n"
        
        # Protein → Pathways
        if pathways:
            description += f"🔹 **{target_gene_symbol}** PARTICIPATES IN:\n"
            for pw in pathways[:3]:
                description += f"   • {pw}\n"
            if len(pathways) > 3:
                description += f"   • ...and {len(pathways)-3} more pathways\n"
            description += "\n"
        
        # Pathways → Disease
        if connected_pathway_count > 0:
            description += f"🔹 **{connected_pathway_count} pathway(s)** LINKED TO **{disease_name}**\n\n"
        else:
            description += f"🔹 Pathways being validated for {disease_name} relevance\n\n"
        
        # Mechanism (from genes.json locus_type or known biology)
        locus_type = gene_info.get('locus_type', '')
        if locus_type:
            description += f"*Mechanism: {target_gene_symbol} is a {locus_type}*"
        
        logger.info(f"✅ Generated evidence chain for {drug_name}")
        return description
        
    except Exception as e:
        logger.error(f"Evidence chain description failed: {e}")
        return f"{drug_name} → {disease_name} evidence chain (details unavailable)"


__all__ = ['generate_evidence_chain_description']
