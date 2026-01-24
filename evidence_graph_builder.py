"""
Evidence Graph Builder - CORRECT Schema
Based on actual table columns provided by user
"""
import psycopg2
import pandas as pd
import logging
import os

logger = logging.getLogger(__name__)

class EvidenceGraphBuilder:
    """Build evidence graphs using CORRECT column names"""
    
    def __init__(self):
        self.conn = None
        self._connect_db()
    
    def _connect_db(self):
        """Connect to database"""
        try:
            self.conn = psycopg2.connect(
                host=os.getenv("DB_HOST", "localhost"),
                database=os.getenv("DB_NAME", "cipherq_repurpose"),
                user=os.getenv("DB_USER", "babburisoumith"),
                password=os.getenv("DB_PASSWORD", "")
            )
            logger.info("✅ EvidenceGraphBuilder connected (200K pathways)")
        except Exception as e:
            logger.error(f"Database connection failed: {e}")
    
    def build_evidence_graph(self, drug_names: list, disease_name: str):
        """Build: Drugs → Proteins → Pathways → Disease"""
        
        if not self.conn:
            return pd.DataFrame(), pd.DataFrame()
        
        try:
            cursor = self.conn.cursor()
            
            nodes = []
            edges = []
            node_ids = set()
            
            # 1. Drug nodes
            for drug in drug_names:
                drug_id = f"DRUG_{drug.replace(' ', '_').upper()}"
                nodes.append({'id': drug_id, 'label': 'Drug', 'name': drug, 'type': 'drug'})
                node_ids.add(drug_id)
            
            # 2. Drugs → Proteins
            # drugs.id → drug_protein_interactions.drug_id
            # proteins.id → drug_protein_interactions.protein_id
            cursor.execute("""
                SELECT d.name, p.gene_symbol, p.name
                FROM drugs d
                JOIN drug_protein_interactions dpi ON d.id = dpi.drug_id
                JOIN proteins p ON p.id = dpi.protein_id
                WHERE d.name = ANY(%s)
                LIMIT 100
            """, (drug_names,))
            
            interactions = cursor.fetchall()
            logger.info(f"✅ {len(interactions)} drug-protein interactions")
            
            protein_db_ids = {}  # gene_symbol → database id
            protein_genes = set()
            
            for drug_name, gene_symbol, protein_name in interactions:
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
                cursor.execute("""
                    SELECT id, gene_symbol FROM proteins WHERE gene_symbol = ANY(%s)
                """, (list(protein_genes),))
                for db_id, gene in cursor.fetchall():
                    protein_db_ids[gene] = db_id
            
            logger.info(f"✅ {len(protein_genes)} proteins")
            
            # 3. Proteins → Pathways  
            # protein_pathway_members: protein_id → pathways.id
            # pathways: has 'name', 'pathway_source' (not pathway_name!)
            if protein_db_ids:
                logger.info(f"🔍 Checking pathways for {len(protein_db_ids)} proteins: {list(protein_db_ids.keys())[:5]}")
                logger.info(f"🔍 Database IDs: {list(protein_db_ids.values())[:5]}")
                
                cursor.execute("""
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
                
                pathway_data = cursor.fetchall()
                
                if len(pathway_data) == 0:
                    logger.warning(f"❌ NO PATHWAYS FOUND for these proteins!")
                    logger.warning(f"   Checked protein IDs: {list(protein_db_ids.values())}")
                    logger.warning(f"   This means protein_pathway_members table is empty or doesn't have these proteins")
                else:
                    logger.info(f"✅ {len(pathway_data)} pathways from 200K!")
                
                pathway_db_ids = set()
                for gene, pw_id, pw_name, pw_source in pathway_data:
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
                # pathway_disease_associations uses disease_id (NOT disease_name!)
                # First get disease_id
                cursor.execute("SELECT id FROM diseases WHERE name ILIKE %s LIMIT 1", (f"%{disease_name}%",))
                disease_row = cursor.fetchone()
                
                if disease_row and pathway_db_ids:
                    disease_db_id = disease_row[0]
                    
                    cursor.execute("""
                        SELECT pathway_id, relevance_score
                        FROM pathway_disease_associations
                        WHERE disease_id = %s
                        AND pathway_id = ANY(%s)
                        LIMIT 20
                    """, (disease_db_id, list(pathway_db_ids)))
                    
                    pw_disease = cursor.fetchall()
                    logger.info(f"✅ {len(pw_disease)} pathway-disease links")
                    
                    disease_id = f"DISEASE_{disease_name.replace(' ', '_').upper()}"
                    for pw_id, score in pw_disease:
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
            
            cursor.close()
            
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
                'diseases': 0
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
            'diseases': diseases_count
        }