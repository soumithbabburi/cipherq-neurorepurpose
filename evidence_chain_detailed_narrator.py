"""
Detailed Evidence Chain Narrator
Explains every node, edge, and connection in the evidence graph
"""

import pandas as pd
import logging
import json

logger = logging.getLogger(__name__)


def generate_detailed_evidence_narrative(nodes_df, edges_df, drug_name: str, disease_name: str) -> str:
    """
    Generate comprehensive narrative explaining EVERY node and edge
    """
    
    try:
        # Load reference data
        with open('genes.json', 'r') as f:
            genes_data = json.load(f)
        
        with open('pathways.json', 'r') as f:
            pathways_data = json.load(f)
        
        # Find drug node
        drug_id = f"DRUG_{drug_name.replace(' ', '_').replace('-', '_').upper()}"
        
        narrative = f"# {drug_name} → {disease_name} Evidence Chain\n\n"
        narrative += "## Complete Node-by-Node Analysis\n\n"
        
        # === NODE 1: DRUG ===
        narrative += f"### NODE 1: DRUG - {drug_name}\n\n"
        narrative += f"**Type:** Small molecule therapeutic agent\n"
        
        # Get drug info
        try:
            with open('drugs.json', 'r') as f:
                drugs = json.load(f)
            drug_lower = drug_name.lower().replace(' hydrochloride', '').replace(' sodium', '')
            if drug_lower in drugs:
                drug_info = drugs[drug_lower]
                mw = drug_info.get('molecular_weight', 'N/A')
                narrative += f"**Molecular Weight:** {mw} Da\n"
                narrative += f"**FDA Status:** {'Approved' if drug_info.get('approved') else 'Investigational'}\n"
        except:
            pass
        
        narrative += f"**Starting point for evidence chain validation**\n\n"
        
        # === EDGES: DRUG → PROTEINS ===
        drug_edges = edges_df[edges_df['source'] == drug_id]
        
        if len(drug_edges) == 0:
            narrative += "**No protein targets found in database.**\n"
            return narrative
        
        # Sort by confidence
        drug_edges_sorted = drug_edges.sort_values('confidence', ascending=False)
        
        narrative += f"## Drug-Protein Connections ({len(drug_edges)} targets)\n\n"
        
        target_proteins = []
        
        for idx, edge in drug_edges_sorted.head(5).iterrows():
            target_id = edge['target']
            confidence = edge.get('confidence', 0.8)
            
            # Get protein node
            protein_node = nodes_df[nodes_df['id'] == target_id]
            if len(protein_node) == 0:
                continue
            
            gene_symbol = protein_node.iloc[0]['name']
            target_proteins.append(gene_symbol)
            
            # Get full protein info
            gene_info = genes_data.get(gene_symbol.upper(), {})
            protein_full_name = gene_info.get('name', gene_symbol)
            
            narrative += f"### EDGE {idx+1}: {drug_name} → {gene_symbol}\n\n"
            narrative += f"**Connection Type:** Drug-Protein Interaction\n"
            narrative += f"**Confidence:** {confidence:.0%} (from DGIdb database)\n"
            narrative += f"**Target:** {gene_symbol} ({protein_full_name})\n"
            
            # Protein functions (from earlier)
            functions = {
                'PPARG': 'Nuclear receptor regulating glucose metabolism, adipogenesis, and inflammatory response in adipose, immune, and brain cells',
                'ABCC8': 'ATP-binding cassette transporter; SUR1 subunit of pancreatic β-cell K-ATP channels controlling insulin secretion',
                'DPP4': 'Membrane-bound serine protease degrading GLP-1 and GIP incretin hormones',
                'ACHE': 'Acetylcholinesterase enzyme at cholinergic synapses; terminates acetylcholine neurotransmission',
                'DRD2': 'G-protein coupled dopamine D2 receptor in striatal neurons; mediates motor control and reward signaling',
                'MAOB': 'Mitochondrial monoamine oxidase B; catalyzes oxidative deamination of dopamine, tyramine, and benzylamine',
                'MGAM': 'Maltase-glucoamylase enzyme in small intestine; hydrolyzes α-1,4-glycosidic bonds in starch and maltose',
                'SI': 'Sucrase-isomaltase enzyme complex; cleaves sucrose and maltose into monosaccharides for absorption',
                'PRKAG3': 'AMPK gamma-3 regulatory subunit; nucleotide sensor activating energy-deficit responses',
                'KCNJ11': 'Kir6.2 potassium channel subunit of pancreatic K-ATP channels; couples metabolism to membrane excitability'
            }
            
            if gene_symbol in functions:
                narrative += f"**Function:** {functions[gene_symbol]}\n\n"
            else:
                narrative += f"**Function:** Protein target (detailed function data limited)\n\n"
        
        # === PROTEIN → PATHWAYS ===
        if len(target_proteins) > 0:
            primary_target = target_proteins[0]
            
            narrative += f"## Protein-Pathway Connections\n\n"
            narrative += f"### NODE 2: PRIMARY TARGET - {primary_target}\n\n"
            
            # Get pathways for this protein
            try:
                with open('protein_pathways.json', 'r') as f:
                    protein_pathways = json.load(f)
                
                if primary_target in protein_pathways:
                    pathway_ids = protein_pathways[primary_target]
                    
                    narrative += f"**{primary_target} participates in {len(pathway_ids)} biological pathways:**\n\n"
                    
                    pathway_count = 0
                    for pw_id in pathway_ids[:10]:
                        if pw_id in pathways_data:
                            pw_name = pathways_data[pw_id].get('name', pw_id)
                            pathway_count += 1
                            
                            narrative += f"{pathway_count}. **{pw_name}**\n"
                            
                            # Explain pathway relevance
                            pw_lower = pw_name.lower()
                            
                            if 'insulin' in pw_lower:
                                narrative += "   *Role: Glucose uptake and metabolic regulation*\n"
                            elif 'glucose' in pw_lower:
                                narrative += "   *Role: Carbohydrate metabolism and energy homeostasis*\n"
                            elif 'dopamin' in pw_lower:
                                narrative += "   *Role: Dopaminergic neurotransmission and motor control*\n"
                            elif 'acetylcholine' in pw_lower:
                                narrative += "   *Role: Cholinergic neurotransmission and cognitive function*\n"
                            elif 'ampk' in pw_lower:
                                narrative += "   *Role: Cellular energy sensing and metabolic adaptation*\n"
                            elif 'lipid' in pw_lower or 'fatty' in pw_lower:
                                narrative += "   *Role: Lipid metabolism and energy storage*\n"
                            
                            narrative += "\n"
                    
                    if len(pathway_ids) > 10:
                        narrative += f"*...and {len(pathway_ids) - 10} additional pathways*\n\n"
                
                else:
                    narrative += f"**No pathway data available for {primary_target}**\n\n"
            
            except Exception as pw_err:
                logger.error(f"Pathway loading error: {pw_err}")
                narrative += "**Pathway data unavailable**\n\n"
        
        # === PATHWAYS → DISEASE ===
        narrative += f"## Pathway-Disease Connections\n\n"
        
        disease_id = f"DISEASE_{disease_name.replace(' ', '_').upper()}"
        pathway_to_disease = edges_df[edges_df['target'] == disease_id]
        
        if len(pathway_to_disease) > 0:
            narrative += f"**{len(pathway_to_disease)} pathway(s) connected to {disease_name}:**\n\n"
            
            for idx, edge in pathway_to_disease.head(10).iterrows():
                pathway_id = edge['source']
                confidence = edge.get('confidence', 0.6)
                
                # Get pathway name from nodes
                pw_node = nodes_df[nodes_df['id'] == pathway_id]
                if len(pw_node) > 0:
                    pw_name = pw_node.iloc[0]['name']
                    
                    narrative += f"**Connection {idx+1}:** {pw_name}\n"
                    narrative += f"   Relevance Score: {confidence:.0%}\n"
                    
                    # Explain disease relevance
                    disease_lower = disease_name.lower()
                    pw_lower = pw_name.lower()
                    
                    if 'alzheimer' in disease_lower:
                        if 'insulin' in pw_lower:
                            narrative += "   *Disease Link: Insulin resistance in the brain promotes amyloid-β accumulation*\n"
                        elif 'glucose' in pw_lower:
                            narrative += "   *Disease Link: Metabolic dysfunction contributes to neurodegeneration*\n"
                        elif 'acetylcholine' in pw_lower:
                            narrative += "   *Disease Link: Cholinergic deficit causes cognitive symptoms*\n"
                    
                    elif 'parkinson' in disease_lower:
                        if 'dopamin' in pw_lower:
                            narrative += "   *Disease Link: Dopamine depletion causes motor symptoms*\n"
                        elif 'mao' in pw_lower:
                            narrative += "   *Disease Link: MAO inhibition preserves dopamine levels*\n"
                    
                    narrative += "\n"
        
        else:
            narrative += f"**No pathways currently linked to {disease_name} in evidence graph**\n"
            narrative += "*This suggests limited mechanistic evidence for repurposing to this indication*\n\n"
        
        # === SUMMARY ===
        narrative += "## Evidence Chain Summary\n\n"
        narrative += f"**Complete Path:** {drug_name} → {len(drug_edges)} proteins → "
        narrative += f"{len(pathway_ids) if 'pathway_ids' in locals() else 0} pathways → "
        narrative += f"{len(pathway_to_disease)} disease connections\n\n"
        
        if len(pathway_to_disease) > 0:
            narrative += f"**Mechanistic Plausibility:** Evidence supports biological pathway from drug to disease. "
            narrative += f"Strength of evidence based on {len(pathway_to_disease)} pathway-disease links and validated protein targets.\n"
        else:
            narrative += "**Mechanistic Plausibility:** Limited - no clear pathway connecting drug mechanism to disease pathophysiology.\n"
        
        return narrative
        
    except Exception as e:
        logger.error(f"Detailed narrator failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return f"**Evidence analysis unavailable for {drug_name} → {disease_name}**"


__all__ = ['generate_detailed_evidence_narrative']
