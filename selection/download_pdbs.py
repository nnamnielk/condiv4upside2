#!/usr/bin/env python3

import os
import sys
import requests
import time
import argparse
from pathlib import Path
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import shutil

def extract_pdb_ids_from_file(info_file):
    """
    Extract unique PDB IDs from a selection info TSV file.
    """
    try:
        df = pd.read_csv(info_file, sep='\t')
        pdb_ids = set()
        
        for pdb_chain in df['PDB_Chain']:
            # Extract PDB ID from entries like "1BII_3|Chain"
            pdb_id = pdb_chain.split('_')[0]
            pdb_ids.add(pdb_id.upper())
        
        return sorted(list(pdb_ids))
    
    except Exception as e:
        print(f"Error reading {info_file}: {e}")
        return []

def download_pdb_file(pdb_id, output_dir, file_format='pdb', compressed=False):
    """
    Download a single PDB file from RCSB PDB.
    
    Args:
        pdb_id: PDB identifier (e.g., '1CRN')
        output_dir: Directory to save the file
        file_format: 'pdb' or 'cif'
        compressed: Whether to download compressed files
    
    Returns:
        tuple: (pdb_id, success, message)
    """
    pdb_id = pdb_id.upper()
    
    # Construct URL based on format
    if file_format == 'pdb':
        if compressed:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb.gz"
            filename = f"{pdb_id}.pdb.gz"
        else:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            filename = f"{pdb_id}.pdb"
    elif file_format == 'cif':
        if compressed:
            url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
            filename = f"{pdb_id}.cif.gz"
        else:
            url = f"https://files.rcsb.org/download/{pdb_id}.cif"
            filename = f"{pdb_id}.cif"
    else:
        return pdb_id, False, f"Unsupported format: {file_format}"
    
    output_path = os.path.join(output_dir, filename)
    
    # Skip if file already exists
    if os.path.exists(output_path):
        return pdb_id, True, f"Already exists: {filename}"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        # Write file
        with open(output_path, 'wb') as f:
            f.write(response.content)
        
        # If compressed and we want to decompress
        if compressed and not compressed:  # This logic needs fixing
            pass  # Keep compressed for now
        
        return pdb_id, True, f"Downloaded: {filename}"
        
    except requests.exceptions.RequestException as e:
        return pdb_id, False, f"Download failed: {str(e)}"
    except Exception as e:
        return pdb_id, False, f"Error: {str(e)}"

def download_pdbs_parallel(pdb_ids, output_dir, file_format='pdb', compressed=False, max_workers=5):
    """
    Download multiple PDB files in parallel.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    successful = []
    failed = []
    
    print(f"Downloading {len(pdb_ids)} PDB files to {output_dir}")
    print(f"Format: {file_format}, Compressed: {compressed}")
    print(f"Using {max_workers} parallel workers")
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all download tasks
        future_to_pdb = {
            executor.submit(download_pdb_file, pdb_id, output_dir, file_format, compressed): pdb_id 
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
            time.sleep(0.1)
    
    return successful, failed

def main():
    parser = argparse.ArgumentParser(description='Download PDB files for selected proteins')
    parser.add_argument('--alpha-info', default='enhanced_alpha_selection_info.tsv',
                       help='Alpha proteins selection info file')
    parser.add_argument('--beta-info', default='enhanced_beta_selection_info.tsv',
                       help='Beta proteins selection info file')
    parser.add_argument('--output-dir', default='pdb_files',
                       help='Output directory for PDB files')
    parser.add_argument('--format', choices=['pdb', 'cif'], default='pdb',
                       help='File format to download')
    parser.add_argument('--compressed', action='store_true',
                       help='Download compressed files (.gz)')
    parser.add_argument('--max-workers', type=int, default=5,
                       help='Maximum number of parallel downloads')
    parser.add_argument('--separate-dirs', action='store_true',
                       help='Create separate directories for alpha and beta proteins')
    parser.add_argument('--alpha-only', action='store_true',
                       help='Download only alpha proteins')
    parser.add_argument('--beta-only', action='store_true',
                       help='Download only beta proteins')
    
    args = parser.parse_args()
    
    # Check if info files exist
    alpha_info_path = args.alpha_info
    beta_info_path = args.beta_info
    
    if not os.path.exists(alpha_info_path) and not args.beta_only:
        print(f"Alpha info file not found: {alpha_info_path}")
        if not args.alpha_only:
            sys.exit(1)
    
    if not os.path.exists(beta_info_path) and not args.alpha_only:
        print(f"Beta info file not found: {beta_info_path}")
        if not args.beta_only:
            sys.exit(1)
    
    all_successful = []
    all_failed = []
    
    # Download alpha proteins
    if not args.beta_only and os.path.exists(alpha_info_path):
        print("\n" + "="*60)
        print("DOWNLOADING ALPHA PROTEINS")
        print("="*60)
        
        alpha_pdb_ids = extract_pdb_ids_from_file(alpha_info_path)
        print(f"Found {len(alpha_pdb_ids)} unique alpha protein PDB IDs")
        
        if args.separate_dirs:
            alpha_output_dir = os.path.join(args.output_dir, 'alpha')
        else:
            alpha_output_dir = args.output_dir
        
        alpha_successful, alpha_failed = download_pdbs_parallel(
            alpha_pdb_ids, alpha_output_dir, args.format, args.compressed, args.max_workers
        )
        
        all_successful.extend(alpha_successful)
        all_failed.extend(alpha_failed)
        
        print(f"\nAlpha proteins: {len(alpha_successful)} successful, {len(alpha_failed)} failed")
    
    # Download beta proteins
    if not args.alpha_only and os.path.exists(beta_info_path):
        print("\n" + "="*60)
        print("DOWNLOADING BETA PROTEINS")
        print("="*60)
        
        beta_pdb_ids = extract_pdb_ids_from_file(beta_info_path)
        print(f"Found {len(beta_pdb_ids)} unique beta protein PDB IDs")
        
        if args.separate_dirs:
            beta_output_dir = os.path.join(args.output_dir, 'beta')
        else:
            beta_output_dir = args.output_dir
        
        beta_successful, beta_failed = download_pdbs_parallel(
            beta_pdb_ids, beta_output_dir, args.format, args.compressed, args.max_workers
        )
        
        all_successful.extend(beta_successful)
        all_failed.extend(beta_failed)
        
        print(f"\nBeta proteins: {len(beta_successful)} successful, {len(beta_failed)} failed")
    
    # Summary
    print("\n" + "="*60)
    print("DOWNLOAD SUMMARY")
    print("="*60)
    print(f"Total successful downloads: {len(all_successful)}")
    print(f"Total failed downloads: {len(all_failed)}")
    
    if all_failed:
        print(f"\nFailed downloads: {', '.join(all_failed)}")
        
        # Write failed IDs to file for retry
        failed_file = os.path.join(args.output_dir, 'failed_downloads.txt')
        with open(failed_file, 'w') as f:
            for pdb_id in all_failed:
                f.write(f"{pdb_id}\n")
        print(f"Failed PDB IDs written to: {failed_file}")
    
    print(f"\nPDB files saved to: {args.output_dir}")

if __name__ == '__main__':
    main()
