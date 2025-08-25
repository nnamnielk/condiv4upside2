#!/usr/bin/env python3

import os
import requests
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

def download_pdb_file(pdb_id, output_dir):
    """
    Download a single PDB file from RCSB PDB.
    """
    pdb_id = pdb_id.upper().strip()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    
    # Skip if file already exists
    if os.path.exists(output_path):
        return pdb_id, True, f"Already exists: {pdb_id}.pdb"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            f.write(response.content)
        
        return pdb_id, True, f"Downloaded: {pdb_id}.pdb"
        
    except requests.exceptions.RequestException as e:
        return pdb_id, False, f"Download failed: {str(e)}"
    except Exception as e:
        return pdb_id, False, f"Error: {str(e)}"

def main():
    missing_file = "missing_pdbs.txt"
    output_dir = "selection/pdbs/alpha"
    max_workers = 3
    
    # Read missing PDB IDs
    if not os.path.exists(missing_file):
        print(f"Missing file {missing_file} not found!")
        return
    
    with open(missing_file, 'r') as f:
        pdb_ids = [line.strip() for line in f if line.strip()]
    
    print(f"Downloading {len(pdb_ids)} missing PDB files to {output_dir}")
    print(f"Using {max_workers} parallel workers")
    
    os.makedirs(output_dir, exist_ok=True)
    
    successful = []
    failed = []
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all download tasks
        future_to_pdb = {
            executor.submit(download_pdb_file, pdb_id, output_dir): pdb_id 
            for pdb_id in pdb_ids
        }
        
        # Process completed downloads
        for i, future in enumerate(as_completed(future_to_pdb), 1):
            pdb_id = future_to_pdb[future]
            try:
                pdb_id, success, message = future.result()
                if success:
                    successful.append(pdb_id)
                    print(f"[{i}/{len(pdb_ids)}] ✓ {message}")
                else:
                    failed.append(pdb_id)
                    print(f"[{i}/{len(pdb_ids)}] ✗ {pdb_id}: {message}")
            except Exception as e:
                failed.append(pdb_id)
                print(f"[{i}/{len(pdb_ids)}] ✗ {pdb_id}: Exception: {e}")
            
            # Small delay to be nice to the server
            time.sleep(0.2)
    
    print(f"\nDownload Summary:")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")
    
    if failed:
        print(f"Failed downloads: {', '.join(failed)}")

if __name__ == '__main__':
    main()
