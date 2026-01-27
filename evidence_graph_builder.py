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

# Load curated interactions
_CURATED_INTERACTIONS = None

def load_curated_interactions():
    """Load curated interactions from JSON"""
    global _CURATED_INTERACTIONS
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

load_curated_interactions()


class EvidenceGraphBuilder:
    """Build evidence graphs using proper BioCypher structure"""
    
    def __init__(self):
        logger.info("EvidenceGraphBuilder initialized")
    
    def build_evidence_graph(self, drug_names: list, disease_name: str):
        """
        Build complete graph: Drugs → Proteins → Genes → Pathways → Disease
        Uses actual database connections - NO HARDCODING
        """
        
        try:
            from database_utils import execute_query
            
            nodes = []
            edges = []
            node_ids = set()
            
            logger.info(f"Building graph for {len(drug_names)} drugs and disease: {disease_name}")
            
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
            
            # 2. GET DRUG → PROTEIN INTERACTIONS (JSON FIRST, THEN DATABASE)
            logger.info("Getting drug-protein interactions...")
            
            interactions = []
            
            # First, try to get from JSON
            for drug_name in drug_names:
                drug_lower = drug_name.lower()
                
                if drug_lower in _CURATED_INTERACTIONS:
                    logger.info(f"✅ {drug_name} found in JSON with {len(_CURATED_INTERACTIONS[drug_lower])} targets")
                    
                    for target in _CURATED_INTERACTIONS[drug_lower]:
                        interactions.append({
                            'drug_name': drug_name,
                            'protein_id': None,  # Not from DB
                            'gene_symbol': target['gene_symbol'],
                            'protein_name': target['protein_name'],
                            'protein_function': '',
                            'confidence_score': target['confidence_score'],
                            'interaction_type': 'target',
                            'binding_affinity': target.get('binding_affinity')
                        })
            
            logger.info(f"Found {len(interactions)} interactions from JSON")
            
            # If no JSON data, try database
            if len(interactions) == 0:
                logger.info("No JSON data found, querying database...")
                
                from database_utils import execute_query
                
                # Check if drugs exist
                drug_check = execute_query("""
                    SELECT name FROM drugs WHERE LOWER(name) = ANY(ARRAY[%s])
                """, ([d.lower() for d in drug_names],))
                
                drugs_found = [d['name'] for d in drug_check]
                
                if not drugs_found:
                    logger.error("Drugs not found in database OR JSON!")
                    return pd.DataFrame([{
                        'id': 'ERROR',
                        'label': 'Error',
                        'name': f'No drugs found. Searched: {drug_names}',
                        'type': 'error'
                    }]), pd.DataFrame()
                
                interactions = execute_query("""
                    SELECT 
                        d.name as drug_name,
                        p.id as protein_id,
                        p.gene_symbol,
                        p.name as protein_name,
                        p.function as protein_function,
                        dpi.confidence_score,
                        dpi.interaction_type,
                        dpi.binding_affinity
                    FROM drugs d
                    JOIN drug_protein_interactions dpi ON d.id = dpi.drug_id
                    JOIN proteins p ON p.id = dpi.protein_id
                    WHERE LOWER(d.name) = ANY(ARRAY[%s])
                    ORDER BY dpi.confidence_score DESC NULLS LAST
                    LIMIT 200
                """, ([d.lower() for d in drugs_found],))
            
            if not interactions:
                logger.warning(f"No interactions found for drugs: {drugs_found}")
                # Try to find ANY interactions for these drugs
                basic_drug_info = execute_query("""
                    SELECT d.name, d.mechanism_of_action, d.drug_class
                    FROM drugs d
                    WHERE LOWER(d.name) = ANY(ARRAY[%s])
                """, ([d.lower() for d in drugs_found],))
                
                logger.info(f"Found basic info for {len(basic_drug_info)} drugs but no interactions")
            
            logger.info(f"Found {len(interactions)} drug-protein interactions")
            
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
                edges.append({
                    'source': drug_id,
                    'target': protein_node_id,
                    'label': interaction_type.upper(),
                    'confidence': float(confidence) if confidence else 0.8,
                    'edge_type': 'drug_protein'
                })
            
            logger.info(f"Added {len(protein_db_ids)} protein/gene nodes")
            
            # 4. GET PROTEIN → PATHWAY CONNECTIONS
            if protein_db_ids:
                logger.info(f"Querying pathways for {len(protein_db_ids)} proteins...")
                
                pathway_rows = execute_query("""
                    SELECT DISTINCT 
                        p.gene_symbol,
                        pw.id as pathway_id,
                        pw.name as pathway_name,
                        pw.pathway_source,
                        pw.pathway_category,
                        ppm.role as protein_role
                    FROM protein_pathway_members ppm
                    JOIN proteins p ON ppm.protein_id = p.id
                    JOIN pathways pw ON ppm.pathway_id = pw.id
                    WHERE ppm.protein_id = ANY(%s)
                    ORDER BY pw.pathway_source, pw.name
                    LIMIT 100
                """, (list(protein_db_ids.values()),))
                
                logger.info(f"Found {len(pathway_rows)} protein-pathway connections")
                
                pathway_db_ids = {}
                
                # 5. ADD PATHWAY NODES
                for row in pathway_rows:
                    gene = row['gene_symbol']
                    pw_id = row['pathway_id']
                    pw_name = row['pathway_name']
                    pw_source = row.get('pathway_source', 'Database')
                    pw_category = row.get('pathway_category', '')
                    protein_role = row.get('protein_role', '')
                    
                    pathway_node_id = f"PATHWAY_{pw_id}"
                    
                    if pathway_node_id not in node_ids:
                        display_name = f"{pw_name}"
                        if pw_category:
                            display_name += f" ({pw_category})"
                        
                        nodes.append({
                            'id': pathway_node_id,
                            'label': 'Pathway',
                            'name': pw_name,
                            'display_name': display_name,
                            'source': pw_source,
                            'category': pw_category,
                            'type': 'pathway',
                            'node_type': 'pathway'
                        })
                        node_ids.add(pathway_node_id)
                        pathway_db_ids[pw_id] = pw_name
                    
                    # Add protein → pathway edge
                    protein_node_id = f"PROTEIN_{gene}"
                    if protein_node_id in node_ids:
                        edges.append({
                            'source': protein_node_id,
                            'target': pathway_node_id,
                            'label': 'PARTICIPATES_IN',
                            'confidence': 0.85,
                            'role': protein_role,
                            'edge_type': 'protein_pathway'
                        })
                
                logger.info(f"Added {len(pathway_db_ids)} pathway nodes")
                
                # 6. CONNECT PATHWAYS → DISEASE
                logger.info(f"Querying disease associations for: {disease_name}")
                
                # First, find the disease in database
                disease_rows = execute_query("""
                    SELECT id, name, disease_category
                    FROM diseases
                    WHERE name ILIKE %s
                    OR disease_category ILIKE %s
                    LIMIT 5
                """, (f"%{disease_name}%", f"%{disease_name}%"))
                
                disease_node_id = f"DISEASE_{disease_name.replace(' ', '_').replace('-', '_').upper()}"
                disease_db_id = None
                disease_category = None
                
                if disease_rows:
                    disease_db_id = disease_rows[0]['id']
                    disease_name_db = disease_rows[0]['name']
                    disease_category = disease_rows[0].get('disease_category')
                    logger.info(f"Found disease in database: {disease_name_db} (ID: {disease_db_id})")
                    
                    # Query pathway-disease associations
                    if pathway_db_ids:
                        pw_disease_rows = execute_query("""
                            SELECT 
                                pda.pathway_id,
                                pda.relevance_score,
                                pda.evidence_source,
                                pw.name as pathway_name
                            FROM pathway_disease_associations pda
                            JOIN pathways pw ON pda.pathway_id = pw.id
                            WHERE pda.disease_id = %s
                            AND pda.pathway_id = ANY(%s)
                            ORDER BY pda.relevance_score DESC NULLS LAST
                            LIMIT 50
                        """, (disease_db_id, list(pathway_db_ids.keys())))
                        
                        logger.info(f"Found {len(pw_disease_rows)} pathway-disease associations")
                        
                        # Add pathway → disease edges
                        for row in pw_disease_rows:
                            pw_id = row['pathway_id']
                            relevance = row.get('relevance_score', 0.75)
                            evidence = row.get('evidence_source', 'Database')
                            
                            pathway_node_id = f"PATHWAY_{pw_id}"
                            if pathway_node_id in node_ids:
                                edges.append({
                                    'source': pathway_node_id,
                                    'target': disease_node_id,
                                    'label': 'IMPLICATED_IN',
                                    'confidence': float(relevance) if relevance else 0.75,
                                    'evidence': evidence,
                                    'edge_type': 'pathway_disease'
                                })
                    else:
                        logger.warning("No pathways to connect to disease")
                else:
                    logger.warning(f"Disease '{disease_name}' not found in database")
            else:
                logger.warning("No proteins found, cannot query pathways")
            
            # 7. ADD DISEASE NODE
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
            
            if not has_pathway_disease_links and protein_genes:
                logger.info("No pathway-disease links found, using direct protein-disease associations")
                
                # Query direct protein-disease associations
                gene_disease_rows = execute_query("""
                    SELECT 
                        p.gene_symbol,
                        gda.association_score,
                        gda.evidence_source
                    FROM gene_disease_associations gda
                    JOIN proteins p ON gda.gene_id = p.id
                    WHERE gda.disease_id = %s
                    AND p.id = ANY(%s)
                    LIMIT 30
                """, (disease_db_id if disease_db_id else 0, list(protein_db_ids.values())))
                
                if gene_disease_rows:
                    logger.info(f"Found {len(gene_disease_rows)} direct gene-disease associations")
                    for row in gene_disease_rows:
                        gene = row['gene_symbol']
                        score = row.get('association_score', 0.7)
                        evidence = row.get('evidence_source', 'Database')
                        
                        protein_node_id = f"PROTEIN_{gene}"
                        if protein_node_id in node_ids:
                            edges.append({
                                'source': protein_node_id,
                                'target': disease_node_id,
                                'label': 'ASSOCIATED_WITH',
                                'confidence': float(score) if score else 0.7,
                                'evidence': evidence,
                                'edge_type': 'protein_disease'
                            })
                else:
                    # Ultimate fallback: connect all proteins to disease with low confidence
                    logger.warning("No direct associations found, using inference")
                    for gene in list(protein_genes.keys())[:10]:
                        protein_node_id = f"PROTEIN_{gene}"
                        edges.append({
                            'source': protein_node_id,
                            'target': disease_node_id,
                            'label': 'POTENTIAL_TARGET',
                            'confidence': 0.5,
                            'evidence': 'Inferred',
                            'edge_type': 'protein_disease_inferred'
                        })
            
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
