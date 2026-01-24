#!/usr/bin/env python3
"""
Migrate drug data from JSON file to PostgreSQL database
"""

import os
import json
import ssl
from pathlib import Path
from urllib.parse import urlparse
import pg8000.native

def migrate_drugs():
    """Migrate drugs from JSON to PostgreSQL using pg8000"""
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
    
    inserted = 0
    skipped = 0
    
    for drug_data in drugs_data:
        try:
            conn.run(
                """
                INSERT INTO drugs (name, drug_class, therapeutic_category, target, mechanism, smiles, source, status)
                VALUES (:name, :drug_class, :therapeutic_category, :target, :mechanism, :smiles, :source, :status)
                ON CONFLICT (name) DO NOTHING
                """,
                name=drug_data.get('name'),
                drug_class=drug_data.get('class'),
                therapeutic_category=drug_data.get('therapeutic_category'),
                target=drug_data.get('target'),
                mechanism=drug_data.get('mechanism'),
                smiles=drug_data.get('smiles', ''),
                source=drug_data.get('source', 'FDA'),
                status=drug_data.get('status', 'Approved')
            )
            inserted += 1
        except Exception as e:
            print(f"Error inserting {drug_data.get('name')}: {e}")
            skipped += 1
    
    result = conn.run("SELECT COUNT(*) FROM drugs")
    count = result[0][0]
    
    print(f"\nMigration complete!")
    print(f"  Inserted: {inserted}")
    print(f"  Skipped (errors): {skipped}")
    print(f"  Total in database: {count}")
    
    conn.close()
    return True

if __name__ == '__main__':
    migrate_drugs()
