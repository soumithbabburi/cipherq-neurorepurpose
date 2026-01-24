#!/usr/bin/env python3
"""
Migrate drug data from JSON file to PostgreSQL database using batch inserts
"""

import os
import json
import ssl
from pathlib import Path
from urllib.parse import urlparse
import pg8000.native

def migrate_drugs_batch():
    """Migrate drugs from JSON to PostgreSQL using batch inserts"""
    db_url = os.environ.get('DATABASE_URL')
    if not db_url:
        print("ERROR: DATABASE_URL not set")
        return False
    
    parsed = urlparse(db_url)
    
    ssl_context = ssl.create_default_context()
    
    try:
        conn = pg8000.native.Connection(
            user=parsed.username,
            password=parsed.password,
            host=parsed.hostname,
            port=parsed.port or 5432,
            database=parsed.path[1:],
            ssl_context=ssl_context
        )
        print("Connected to PostgreSQL database")
    except Exception as e:
        print(f"Connection error: {e}")
        return False
    
    json_path = Path(__file__).parent.parent / 'data' / 'drugs_40k.json'
    if not json_path.exists():
        print(f"ERROR: JSON file not found at {json_path}")
        return False
    
    with open(json_path, 'r') as f:
        drugs_data = json.load(f)
    
    print(f"Loading {len(drugs_data)} drugs from JSON...")
    
    conn.run("DELETE FROM drugs")
    print("Cleared existing drugs table")
    
    batch_size = 100
    total_inserted = 0
    
    for i in range(0, len(drugs_data), batch_size):
        batch = drugs_data[i:i+batch_size]
        
        values_list = []
        params = {}
        
        for j, drug in enumerate(batch):
            idx = i + j
            values_list.append(f"(:name_{idx}, :drug_class_{idx}, :therapeutic_category_{idx}, :target_{idx}, :mechanism_{idx}, :smiles_{idx}, :source_{idx}, :status_{idx})")
            params[f'name_{idx}'] = drug.get('name')
            params[f'drug_class_{idx}'] = drug.get('class')
            params[f'therapeutic_category_{idx}'] = drug.get('therapeutic_category')
            params[f'target_{idx}'] = drug.get('target')
            params[f'mechanism_{idx}'] = drug.get('mechanism')
            params[f'smiles_{idx}'] = drug.get('smiles', '')
            params[f'source_{idx}'] = drug.get('source', 'FDA')
            params[f'status_{idx}'] = drug.get('status', 'Approved')
        
        values_sql = ", ".join(values_list)
        sql = f"""
            INSERT INTO drugs (name, drug_class, therapeutic_category, target, mechanism, smiles, source, status)
            VALUES {values_sql}
            ON CONFLICT (name) DO NOTHING
        """
        
        try:
            conn.run(sql, **params)
            total_inserted += len(batch)
            print(f"  Inserted batch {i//batch_size + 1}: {len(batch)} drugs (total: {total_inserted})")
        except Exception as e:
            print(f"Error inserting batch: {e}")
    
    result = conn.run("SELECT COUNT(*) FROM drugs")
    count = result[0][0]
    
    result_categories = conn.run("""
        SELECT therapeutic_category, COUNT(*) as cnt 
        FROM drugs 
        GROUP BY therapeutic_category 
        ORDER BY cnt DESC
    """)
    
    print(f"\nMigration complete!")
    print(f"  Total in database: {count}")
    print(f"\nDrugs by category:")
    for row in result_categories:
        print(f"  {row[0]}: {row[1]}")
    
    conn.close()
    return True

if __name__ == '__main__':
    migrate_drugs_batch()
