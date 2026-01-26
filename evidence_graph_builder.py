"""
Evidence Graph Builder - Uses Centralized Database Utils
No more "connection already closed" errors!
"""
import pandas as pd
import logging
import os
import sys

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

logger = logging.getLogger(__name__)

class EvidenceGraphBuilder:
    """Build evidence graphs using centralized database connection"""
    
    def __init__(self):
        # Don't create connection - use database_utils!
        logger.info("✅ EvidenceGraphBuilder initialized (using database_utils)")
    
    def build_evidence_graph(self, drug_names: list, disease_name: str):
        """Build: Drugs → Proteins → Pathways → Disease"""
        
        try:
            from database_utils import execute_query
            
            nodes = []
            edges = []
            node_ids = set()
            
            # 1. Drug nodes
            for drug in drug_names:
                drug_id = f"DRUG_{drug.replace(' ', '_').upper()}"
                nodes.append({'id': drug_id, 'label': 'Drug', 'name': drug, 'type': 'drug'})
                node_ids.add(drug_id)
            
            # 2. Drugs → Proteins (uses database_utils)
            interactions = execute_query("""
                SELECT d.name, p.gene_symbol, p.name
                FROM drugs d
                JOIN drug_protein_interactions dpi ON d.id = dpi.drug_id
                JOIN proteins p ON p.id = dpi.protein_id
                WHERE d.name = ANY(%s)
                LIMIT 100
            """, (drug_names,))
            
            logger.info(f"✅ {len(interactions)} drug-protein interactions")
            
            protein_db_ids = {}
            protein_genes = set()
            
            for row in interactions:
                drug_name = row['name']
                gene_symbol = row['gene_symbol']
                protein_name = row.get('name')
                
                protein_id = f"PROTEIN_{gene_symbol}"
                if protein_id not in node_ids:
                    nodes.append({
                        'id': protein_id,
                        'label': 'Protein',
                        'name': protein_name or gene_symbol,
                        'type': 'protein'
                    })
                    node_ids.add(protein_id)
                    protein_genes.add(gene_symbol)
                
                drug_id = f"DRUG_{drug_name.replace(' ', '_').upper()}"
                edges.append({'source': drug_id, 'target': protein_id, 'label': 'TARGETS', 'confidence': 0.9})
            
            # Get database IDs for proteins
            if protein_genes:
                protein_rows = execute_query("""
                    SELECT id, gene_symbol FROM proteins WHERE gene_symbol = ANY(%s)
                """, (list(protein_genes),))
                
                for row in protein_rows:
                    protein_db_ids[row['gene_symbol']] = row['id']
            
            logger.info(f"✅ {len(protein_genes)} proteins")
            
            # 3. Proteins → Pathways
            if protein_db_ids:
                logger.info(f"🔍 Checking pathways for {len(protein_db_ids)} proteins")
                
                pathway_rows = execute_query("""
                    SELECT DISTINCT 
                        p.gene_symbol,
                        pw.id,
                        pw.name,
                        pw.pathway_source
                    FROM protein_pathway_members ppm
                    JOIN proteins p ON ppm.protein_id = p.id
                    JOIN pathways pw ON ppm.pathway_id = pw.id
                    WHERE ppm.protein_id = ANY(%s)
                    LIMIT 30
                """, (list(protein_db_ids.values()),))
                
                if len(pathway_rows) == 0:
                    logger.warning(f"❌ NO PATHWAYS (protein_pathway_members is empty)")
                else:
                    logger.info(f"✅ {len(pathway_rows)} pathways from 200K!")
                
                pathway_db_ids = set()
                for row in pathway_rows:
                    gene = row['gene_symbol']
                    pw_id = row['id']
                    pw_name = row['name']
                    pw_source = row.get('pathway_source')
                    
                    pathway_node_id = f"PATHWAY_{pw_id}"
                    
                    if pathway_node_id not in node_ids:
                        nodes.append({
                            'id': pathway_node_id,
                            'label': 'Pathway',
                            'name': f"{pw_name} ({pw_source})" if pw_source else pw_name,
                            'type': 'pathway'
                        })
                        node_ids.add(pathway_node_id)
                        pathway_db_ids.add(pw_id)
                    
                    protein_id = f"PROTEIN_{gene}"
                    edges.append({'source': protein_id, 'target': pathway_node_id, 'label': 'PARTICIPATES_IN', 'confidence': 0.85})
                
                # 4. Pathways → Disease
                disease_rows = execute_query("SELECT id FROM diseases WHERE name ILIKE %s LIMIT 1", (f"%{disease_name}%",))
                
                if disease_rows and pathway_db_ids:
                    disease_db_id = disease_rows[0]['id']
                    
                    pw_disease_rows = execute_query("""
                        SELECT pathway_id, relevance_score
                        FROM pathway_disease_associations
                        WHERE disease_id = %s
                        AND pathway_id = ANY(%s)
                        LIMIT 20
                    """, (disease_db_id, list(pathway_db_ids)))
                    
                    logger.info(f"✅ {len(pw_disease_rows)} pathway-disease links")
                    
                    disease_id = f"DISEASE_{disease_name.replace(' ', '_').upper()}"
                    for row in pw_disease_rows:
                        pw_id = row['pathway_id']
                        score = row.get('relevance_score')
                        pathway_node_id = f"PATHWAY_{pw_id}"
                        if pathway_node_id in node_ids:
                            edges.append({'source': pathway_node_id, 'target': disease_id, 'label': 'IMPLICATED_IN', 'confidence': score or 0.8})
            
            # 5. Disease node
            disease_id = f"DISEASE_{disease_name.replace(' ', '_').upper()}"
            nodes.append({'id': disease_id, 'label': 'Disease', 'name': disease_name, 'type': 'disease'})
            
            # Fallback: connect proteins directly if no pathway links
            if not any(e.get('target') == disease_id for e in edges):
                for gene in protein_genes:
                    protein_id = f"PROTEIN_{gene}"
                    edges.append({'source': protein_id, 'target': disease_id, 'label': 'ASSOCIATED', 'confidence': 0.7})
            
            # Convert to DataFrames
            nodes_df = pd.DataFrame(nodes)
            edges_df = pd.DataFrame(edges)
            
            logger.info(f"✅ Graph: {len(nodes_df)} nodes, {len(edges_df)} edges")
            logger.info(f"   Drugs: {sum(1 for n in nodes if n['type']=='drug')}, Proteins: {sum(1 for n in nodes if n['type']=='protein')}, Pathways: {sum(1 for n in nodes if n['type']=='pathway')}, Disease: {sum(1 for n in nodes if n['type']=='disease')}")
            
            return nodes_df, edges_df
            
        except Exception as e:
            logger.error(f"Graph building failed: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return pd.DataFrame(), pd.DataFrame()
    
    def get_summary_metrics(self, nodes_df, edges_df):
        """Get summary metrics for the graph"""
        
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
            'pathways': pathways_count,
            'diseases': diseases_count,
            'target_count': proteins_count
        }
