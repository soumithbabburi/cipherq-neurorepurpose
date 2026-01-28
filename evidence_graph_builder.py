"""
BioCypher Evidence Graph Builder - FIXED VERSION
Properly connects: Drugs → Proteins/Genes → Pathways → Disease
Checks drug_interactions.json FIRST, then database
"""
import pandas as pd
import logging
import os
import sys
import json

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

logger = logging.getLogger(__name__)

# Load curated interactions and pathways and genes
_CURATED_INTERACTIONS = None
_PATHWAYS = None
_PROTEIN_PATHWAYS = None
_GENES = None

def load_curated_data():
    """Load curated interactions, pathway data, and official genes from JSON"""
    global _CURATED_INTERACTIONS, _PATHWAYS, _PROTEIN_PATHWAYS, _GENES
    
    if _CURATED_INTERACTIONS is not None:
        return
    
    try:
        if os.path.exists('drug_interactions.json'):
            with open('drug_interactions.json', 'r') as f:
                _CURATED_INTERACTIONS = json.load(f)
            logger.info(f"✅ Loaded curated interactions for {len(_CURATED_INTERACTIONS)} drugs")
        else:
            _CURATED_INTERACTIONS = {}
    except:
        _CURATED_INTERACTIONS = {}
    
    try:
        if os.path.exists('pathways.json'):
            with open('pathways.json', 'r') as f:
                _PATHWAYS = json.load(f)
            logger.info(f"✅ Loaded {len(_PATHWAYS)} pathways")
        else:
            _PATHWAYS = {}
    except:
        _PATHWAYS = {}
    
    try:
        if os.path.exists('protein_pathways.json'):
            with open('protein_pathways.json', 'r') as f:
                _PROTEIN_PATHWAYS = json.load(f)
            logger.info(f"✅ Loaded protein-pathway mappings")
        else:
            _PROTEIN_PATHWAYS = {}
    except:
        _PROTEIN_PATHWAYS = {}
    
    try:
        if os.path.exists('genes.json'):
            with open('genes.json', 'r') as f:
                _GENES = json.load(f)
            logger.info(f"✅ Loaded {len(_GENES)} official HGNC genes")
        else:
            _GENES = {}
    except:
        _GENES = {}

load_curated_data()


class EvidenceGraphBuilder:
    """Build evidence graphs using proper BioCypher structure"""
    
    def __init__(self):
        logger.info("EvidenceGraphBuilder initialized")
    
    def build_evidence_graph(self, drug_names: list, disease_name: str):
        """
        Build complete graph: Drugs → Proteins → Pathways → Disease
        Uses JSON files ONLY - NO DATABASE
        """
        
        try:
            # Load JSON data
            load_curated_data()
            
            nodes = []
            edges = []
            node_ids = set()
            
            logger.info(f"Building graph for {len(drug_names)} drugs and disease: {disease_name}")
            logger.info(f"Drug names to search: {drug_names}")
            
            # 1. ADD DRUG NODES
            for drug in drug_names:
                drug_id = f"DRUG_{drug.replace(' ', '_').replace('-', '_').upper()}"
                nodes.append({
                    'id': drug_id,
                    'label': 'Drug',
                    'name': drug,
                    'type': 'drug',
                    'node_type': 'drug'
                })
                node_ids.add(drug_id)
            
            # 2. GET DRUG → PROTEIN INTERACTIONS (JSON ONLY)
            logger.info("Getting drug-protein interactions from JSON...")
            
            interactions = []
            
            # Try to get from JSON with FLEXIBLE name matching
            for drug_name in drug_names:
                drug_lower = drug_name.lower()
                
                # Try exact match first
                if drug_lower in _CURATED_INTERACTIONS:
                    targets_found = _CURATED_INTERACTIONS[drug_lower]
                    logger.info(f"✅ {drug_name} EXACT match in JSON")
                else:
                    # Try partial match: "Pioglitazone" matches "pioglitazone hydrochloride"
                    targets_found = None
                    for json_key in _CURATED_INTERACTIONS.keys():
                        if drug_lower in json_key or json_key in drug_lower:
                            targets_found = _CURATED_INTERACTIONS[json_key]
                            logger.info(f"✅ {drug_name} PARTIAL match: '{json_key}'")
                            break
                
                if not targets_found:
                    logger.warning(f"⚠️ {drug_name} NOT FOUND in JSON!")
                    logger.warning(f"   Available drugs (sample): {list(_CURATED_INTERACTIONS.keys())[:10]}")
                    continue
                
                all_targets = targets_found
                logger.info(f"   Found {len(all_targets)} targets for {drug_name}")
                
                # SORT by confidence and take TOP 5 only!
                sorted_targets = sorted(all_targets, key=lambda x: x['confidence_score'], reverse=True)
                top_targets = sorted_targets[:5]
                
                logger.info(f"   Using top {len(top_targets)} targets (confidence: {top_targets[0]['confidence_score']:.2f} to {top_targets[-1]['confidence_score']:.2f})")
                
                for target in top_targets:
                        gene_symbol = target['gene_symbol']
                        
                        # Get official protein name from genes.json if available
                        protein_name = target.get('protein_name', gene_symbol)
                        if _GENES and gene_symbol in _GENES:
                            official_name = _GENES[gene_symbol].get('name', protein_name)
                            if official_name and official_name != gene_symbol:
                                protein_name = official_name
                        
                        interactions.append({
                            'drug_name': drug_name,
                            'protein_id': None,  # Not from DB
                            'gene_symbol': gene_symbol,
                            'protein_name': protein_name,
                            'protein_function': '',
                            'confidence_score': target['confidence_score'],
                            'interaction_type': target.get('interaction_type', 'target'),
                            'binding_affinity': target.get('binding_affinity')
                        })
            
            logger.info(f"Found {len(interactions)} interactions from JSON")
            
            # If no JSON data, return error explaining the issue
            if len(interactions) == 0:
                logger.error(f"No interactions found in JSON for drugs: {drug_names}")
                logger.error(f"JSON file has these drugs: {list(_CURATED_INTERACTIONS.keys())[:20]}")
                logger.error(f"Searched for (lowercase): {[d.lower() for d in drug_names]}")
                
                return pd.DataFrame([{
                    'id': 'ERROR',
                    'label': 'Error',
                    'name': f'No drug-protein data found in JSON for: {", ".join(drug_names)}',
                    'type': 'error'
                }]), pd.DataFrame()
            
            if not interactions:
                logger.warning(f"No interactions found in JSON")
                return pd.DataFrame([{
                    'id': 'ERROR',
                    'label': 'Error',
                    'name': f'No interactions in JSON',
                    'type': 'error'
                }]), pd.DataFrame()
            
            logger.info(f"Found {len(interactions)} drug-protein interactions from JSON")
            
            protein_db_ids = {}
            protein_genes = {}
            
            # 3. ADD PROTEIN/GENE NODES
            for row in interactions:
                drug_name = row['drug_name']
                protein_id_db = row['protein_id']
                gene_symbol = row['gene_symbol']
                protein_name = row.get('protein_name')
                protein_function = row.get('protein_function', '')
                confidence = row.get('confidence_score', 0.8)
                interaction_type = row.get('interaction_type', 'targets')
                
                # Create protein/gene node
                protein_node_id = f"PROTEIN_{gene_symbol}"
                if protein_node_id not in node_ids:
                    nodes.append({
                        'id': protein_node_id,
                        'label': 'Protein/Gene',
                        'name': f"{gene_symbol}",
                        'full_name': protein_name or gene_symbol,
                        'function': protein_function[:100] if protein_function else 'Protein target',
                        'type': 'protein',
                        'node_type': 'protein'
                    })
                    node_ids.add(protein_node_id)
                    protein_db_ids[gene_symbol] = protein_id_db
                    protein_genes[gene_symbol] = {
                        'name': protein_name,
                        'function': protein_function
                    }
                
                # Add drug → protein edge
                drug_id = f"DRUG_{drug_name.replace(' ', '_').replace('-', '_').upper()}"
                
                # DEBUG: Make sure drug_id exists in node_ids
                if drug_id not in node_ids:
                    logger.warning(f"Drug ID {drug_id} not in node_ids! Adding it now...")
                    # Add missing drug node
                    nodes.append({
                        'id': drug_id,
                        'label': 'Drug',
                        'name': drug_name,
                        'type': 'drug',
                        'node_type': 'drug'
                    })
                    node_ids.add(drug_id)
                
                edges.append({
                    'source': drug_id,
                    'target': protein_node_id,
                    'label': 'TARGETS',
                    'confidence': float(confidence) if confidence else 0.8,
                    'edge_type': 'drug_protein'
                })
                
                logger.debug(f"Added edge: {drug_id} → {protein_node_id}")
            
            logger.info(f"Added {len(protein_db_ids)} protein/gene nodes")
            logger.info(f"Added {sum(1 for e in edges if e.get('edge_type')=='drug_protein')} drug-protein edges")
            
            # 4. GET PROTEIN → PATHWAY CONNECTIONS FROM JSON
            logger.info(f"Getting pathways for {len(protein_genes)} proteins from JSON...")
            
            pathway_nodes_added = {}
            
            if _PROTEIN_PATHWAYS and _PATHWAYS:
                for gene_symbol in protein_genes.keys():
                    # Get pathways for this protein from JSON
                    pathway_ids = _PROTEIN_PATHWAYS.get(gene_symbol, [])
                    
                    for pathway_id in pathway_ids:
                        pathway_info = _PATHWAYS.get(pathway_id)
                        if not pathway_info:
                            continue
                        
                        pathway_node_id = f"PATHWAY_{pathway_id}"
                        
                        # Add pathway node if not already added
                        if pathway_node_id not in node_ids:
                            nodes.append({
                                'id': pathway_node_id,
                                'label': 'Pathway',
                                'name': pathway_info['name'],
                                'display_name': f"{pathway_info['name']} ({pathway_info.get('category', '')})",
                                'source': pathway_info.get('source', 'KEGG'),
                                'category': pathway_info.get('category', ''),
                                'type': 'pathway',
                                'node_type': 'pathway'
                            })
                            node_ids.add(pathway_node_id)
                            pathway_nodes_added[pathway_id] = pathway_info
                        
                        # Add protein → pathway edge
                        protein_node_id = f"PROTEIN_{gene_symbol}"
                        if protein_node_id in node_ids:
                            edges.append({
                                'source': protein_node_id,
                                'target': pathway_node_id,
                                'label': 'PARTICIPATES_IN',
                                'confidence': 0.85,
                                'edge_type': 'protein_pathway'
                            })
                
                logger.info(f"Added {len(pathway_nodes_added)} pathway nodes from JSON")
                
                # 5. ADD DISEASE NODE (if not already added) AND CONNECT PATHWAYS
                disease_node_id = f"DISEASE_{disease_name.replace(' ', '_').replace('-', '_').upper()}"
                
                if disease_node_id not in node_ids:
                    nodes.append({
                        'id': disease_node_id,
                        'label': 'Disease',
                        'name': disease_name,
                        'type': 'disease',
                        'node_type': 'disease'
                    })
                    node_ids.add(disease_node_id)
                
                # Connect pathways to disease using JSON pathway data
                pathway_disease_connections = 0
                
                for pathway_id, pathway_info in pathway_nodes_added.items():
                    pathway_diseases = pathway_info.get('diseases', [])
                    pathway_name = pathway_info.get('name', '')
                    relevance_score = pathway_info.get('relevance_score', 0.7)
                    
                    # Check if this pathway is relevant to the target disease
                    disease_match = False
                    
                    # Match against pathway's disease list
                    for disease_item in pathway_diseases:
                        if disease_name.lower() in disease_item.lower() or disease_item.lower() in disease_name.lower():
                            disease_match = True
                            break
                    
                    # Also check pathway NAME for disease keywords
                    if not disease_match:
                        disease_keywords = disease_name.lower().split()
                        pathway_lower = pathway_name.lower()
                        if any(keyword in pathway_lower for keyword in disease_keywords if len(keyword) > 3):
                            disease_match = True
                    
                    if disease_match:
                        pathway_node_id = f"PATHWAY_{pathway_id}"
                        edges.append({
                            'source': pathway_node_id,
                            'target': disease_node_id,
                            'label': 'ASSOCIATED_WITH',
                            'confidence': relevance_score,
                            'edge_type': 'pathway_disease'
                        })
                        pathway_disease_connections += 1
                
                logger.info(f"Connected {pathway_disease_connections} pathways to {disease_name}")
            else:
                logger.warning("No pathway data loaded from JSON")
            
            # NO DATABASE FALLBACK - JSON ONLY!
            # All pathway and disease data comes from JSON files
            
            # 7. ADD DISEASE NODE (if not already added)
            disease_node_id = f"DISEASE_{disease_name.replace(' ', '_').replace('-', '_').upper()}"
            
            if disease_node_id not in node_ids:
                disease_display = disease_name
                if disease_category:
                    disease_display += f" ({disease_category})"
                
                nodes.append({
                    'id': disease_node_id,
                    'label': 'Disease',
                    'name': disease_name,
                    'display_name': disease_display,
                    'category': disease_category,
                    'type': 'disease',
                    'node_type': 'disease'
                })
                node_ids.add(disease_node_id)
            
            # 8. FALLBACK: If no pathway-disease connections, connect proteins directly
            has_pathway_disease_links = any(
                e.get('edge_type') == 'pathway_disease' for e in edges
            )
            
            # NO DATABASE FALLBACK - all connections from JSON!
            # Pathway-disease links already added from pathways.json
            
            # Also need to define disease_category variable
            disease_category = None
            
            # 9. CONVERT TO DATAFRAMES
            nodes_df = pd.DataFrame(nodes)
            edges_df = pd.DataFrame(edges)
            
            # Log final statistics
            drug_count = sum(1 for n in nodes if n.get('type') == 'drug')
            protein_count = sum(1 for n in nodes if n.get('type') == 'protein')
            pathway_count = sum(1 for n in nodes if n.get('type') == 'pathway')
            disease_count = sum(1 for n in nodes if n.get('type') == 'disease')
            
            logger.info(f"GRAPH BUILT: {len(nodes)} nodes, {len(edges)} edges")
            logger.info(f"  Drugs: {drug_count}, Proteins: {protein_count}, Pathways: {pathway_count}, Disease: {disease_count}")
            logger.info(f"  Edges: Drug→Protein: {sum(1 for e in edges if e.get('edge_type')=='drug_protein')}, "
                       f"Protein→Pathway: {sum(1 for e in edges if e.get('edge_type')=='protein_pathway')}, "
                       f"Pathway→Disease: {sum(1 for e in edges if e.get('edge_type')=='pathway_disease')}")
            
            return nodes_df, edges_df
            
        except Exception as e:
            logger.error(f"Graph building failed: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return pd.DataFrame(), pd.DataFrame()
    
    def get_summary_metrics(self, nodes_df, edges_df):
        """Get comprehensive graph metrics"""
        
        if nodes_df.empty or edges_df.empty:
            return {
                'total_nodes': 0,
                'total_edges': 0,
                'drug_count': 0,
                'protein_count': 0,
                'pathway_count': 0,
                'disease_count': 0,
                'drugs': 0,
                'proteins': 0,
                'genes': 0,
                'pathways': 0,
                'diseases': 0,
                'target_count': 0
            }
        
        drugs_count = len(nodes_df[nodes_df['type'] == 'drug']) if 'type' in nodes_df.columns else 0
        proteins_count = len(nodes_df[nodes_df['type'] == 'protein']) if 'type' in nodes_df.columns else 0
        pathways_count = len(nodes_df[nodes_df['type'] == 'pathway']) if 'type' in nodes_df.columns else 0
        diseases_count = len(nodes_df[nodes_df['type'] == 'disease']) if 'type' in nodes_df.columns else 0
        
        return {
            'total_nodes': len(nodes_df),
            'total_edges': len(edges_df),
            'drug_count': drugs_count,
            'protein_count': proteins_count,
            'pathway_count': pathways_count,
            'disease_count': diseases_count,
            'drugs': drugs_count,
            'proteins': proteins_count,
            'genes': proteins_count,
            'pathways': pathways_count,
            'diseases': diseases_count,
            'target_count': proteins_count
        }
    
    def get_pathway_details(self, nodes_df, edges_df):
        """Extract detailed pathway information for display"""
        
        if nodes_df.empty or edges_df.empty:
            return []
        
        # Get all pathway nodes
        pathway_nodes = nodes_df[nodes_df['type'] == 'pathway'].to_dict('records') if 'type' in nodes_df.columns else []
        
        pathway_details = []
        for pathway in pathway_nodes:
            pathway_id = pathway['id']
            
            # Find proteins in this pathway
            pathway_edges = edges_df[
                (edges_df['target'] == pathway_id) & 
                (edges_df.get('edge_type', '') == 'protein_pathway')
            ] if 'target' in edges_df.columns else pd.DataFrame()
            
            proteins_in_pathway = []
            if not pathway_edges.empty:
                for _, edge in pathway_edges.iterrows():
                    protein_id = edge['source']
                    protein_node = nodes_df[nodes_df['id'] == protein_id]
                    if not protein_node.empty:
                        proteins_in_pathway.append(protein_node.iloc[0]['name'])
            
            # Check if connected to disease
            disease_edges = edges_df[
                (edges_df['source'] == pathway_id) &
                (edges_df.get('edge_type', '').str.contains('disease', na=False))
            ] if 'source' in edges_df.columns else pd.DataFrame()
            
            is_disease_relevant = not disease_edges.empty
            relevance_score = disease_edges.iloc[0]['confidence'] if not disease_edges.empty and 'confidence' in disease_edges.columns else 0.0
            
            pathway_details.append({
                'name': pathway.get('name', 'Unknown'),
                'display_name': pathway.get('display_name', pathway.get('name', 'Unknown')),
                'source': pathway.get('source', 'Database'),
                'category': pathway.get('category', ''),
                'proteins': proteins_in_pathway,
                'protein_count': len(proteins_in_pathway),
                'disease_relevant': is_disease_relevant,
                'relevance_score': relevance_score
            })
        
        # Sort by relevance
        pathway_details.sort(key=lambda x: (x['disease_relevant'], x['relevance_score']), reverse=True)
        
        return pathway_details
