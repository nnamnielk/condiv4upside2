#!/usr/bin/env python3

import os
import requests
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

def download_pdb_file(pdb_id, output_dir):
    """Download a single PDB file from RCSB PDB."""
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
    pdb_list_file = "final_huge_pdb_list.txt"
    output_dir = "selection/pdbs/alpha"
    max_workers = 3
    
    # Read PDB IDs
    with open(pdb_list_file, 'r') as f:
        pdb_ids = [line.strip() for line in f if line.strip()]
    
    print(f"Downloading {len(pdb_ids)} PDB files for optimal training set")
    print(f"Output directory: {output_dir}")
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
                    if i % 10 == 0 or i <= 20:
                        print(f"[{i}/{len(pdb_ids)}] ✓ {message}")
                else:
                    failed.append(pdb_id)
                    print(f"[{i}/{len(pdb_ids)}] ✗ {pdb_id}: {message}")
            except Exception as e:
                failed.append(pdb_id)
                print(f"[{i}/{len(pdb_ids)}] ✗ {pdb_id}: Exception: {e}")
            
            # Small delay to be nice to the server
            time.sleep(0.1)
    
    print(f"\nFinal Download Summary:")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")
    print(f"Success rate: {len(successful)/len(pdb_ids)*100:.1f}%")
    
    if failed:
        print(f"Failed downloads: {', '.join(failed[:10])}{'...' if len(failed) > 10 else ''}")

if __name__ == '__main__':
    main()
