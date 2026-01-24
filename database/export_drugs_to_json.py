#!/usr/bin/env python3
"""
Export drugs from PostgreSQL database to JSON for local deployment.
Creates portable data files that work without database connection.
"""

import os
import sys
import json
import csv
from pathlib import Path
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

def get_db_connection():
    """Get database connection"""
    try:
        import pg8000
        db_url = os.environ.get('DATABASE_URL', '')
        if not db_url:
            print("DATABASE_URL not set")
            return None
        
        # Parse connection string
        if db_url.startswith('postgresql://'):
            db_url = db_url.replace('postgresql://', '')
        
        user_pass, host_db = db_url.split('@')
        user, password = user_pass.split(':')
        host_port, database = host_db.split('/')
        
        if ':' in host_port:
            host, port = host_port.split(':')
            port = int(port)
        else:
            host = host_port
            port = 5432
        
        conn = pg8000.connect(
            user=user,
            password=password,
            host=host,
            port=port,
            database=database,
            ssl_context=True
        )
        return conn
    except Exception as e:
        print(f"Database connection failed: {e}")
        return None

def export_drugs():
    """Export all drugs to JSON and CSV"""
    conn = get_db_connection()
    if not conn:
        print("Cannot connect to database. Using existing JSON files.")
        return False
    
    try:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT name, drug_class, therapeutic_category, target, mechanism, smiles, source, status
            FROM drugs
            ORDER BY name
        """)
        
        rows = cursor.fetchall()
        
        drugs = []
        for row in rows:
            drugs.append({
                'name': row[0],
                'class': row[1],
                'therapeutic_category': row[2],
                'target': row[3],
                'mechanism': row[4],
                'smiles': row[5] or '',
                'source': row[6],
                'status': row[7]
            })
        
        conn.close()
        
        # Create export directory
        export_dir = Path(__file__).parent.parent / 'data_exports'
        export_dir.mkdir(exist_ok=True)
        
        # Export to JSON
        json_path = export_dir / 'drugs_full_export.json'
        with open(json_path, 'w') as f:
            json.dump(drugs, f, indent=2)
        print(f"Exported {len(drugs)} drugs to {json_path}")
        
        # Export to CSV
        csv_path = export_dir / 'drugs_full_export.csv'
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['name', 'class', 'therapeutic_category', 'target', 'mechanism', 'smiles', 'source', 'status'])
            writer.writeheader()
            writer.writerows(drugs)
        print(f"Exported {len(drugs)} drugs to {csv_path}")
        
        # Also update data/drugs_40k.json for loader compatibility
        data_dir = Path(__file__).parent.parent / 'data'
        drugs_40k_path = data_dir / 'drugs_40k.json'
        with open(drugs_40k_path, 'w') as f:
            json.dump(drugs, f, indent=2)
        print(f"Updated {drugs_40k_path}")
        
        # Create metadata file
        metadata = {
            'export_date': datetime.now().isoformat(),
            'total_drugs': len(drugs),
            'categories': {},
            'sources': {}
        }
        
        for drug in drugs:
            cat = drug.get('therapeutic_category', 'unknown')
            src = drug.get('source', 'unknown')
            metadata['categories'][cat] = metadata['categories'].get(cat, 0) + 1
            metadata['sources'][src] = metadata['sources'].get(src, 0) + 1
        
        meta_path = export_dir / 'export_metadata.json'
        with open(meta_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        print(f"Created metadata at {meta_path}")
        
        return True
        
    except Exception as e:
        print(f"Export error: {e}")
        if conn:
            conn.close()
        return False

if __name__ == '__main__':
    print("CipherQ Drug Database Export")
    print("=" * 40)
    success = export_drugs()
    if success:
        print("\nExport completed successfully!")
        print("Files created in data_exports/:")
        print("  - drugs_full_export.json")
        print("  - drugs_full_export.csv")
        print("  - export_metadata.json")
    else:
        print("\nExport failed or database unavailable.")
        print("Use existing JSON files in data/ folder.")
