#!/usr/bin/env python3

import os
import shutil
import pandas as pd
import subprocess
import sys

def get_smallest_alpha_proteins(selection_file, n=15):
    """
    Get the n smallest alpha proteins from the selection info file.
    """
    df = pd.read_csv(selection_file, sep='\t')
    
    # Sort by length and get the smallest n proteins
    df_sorted = df.sort_values('Length').head(n)
    
    # Extract PDB IDs
    pdb_ids = []
    for pdb_chain in df_sorted['PDB_Chain']:
        pdb_id = pdb_chain.split('_')[0]
        pdb_ids.append(pdb_id)
    
    return pdb_ids, df_sorted

def copy_pdb_files(pdb_ids, source_dir, dest_dir):
    """
    Copy PDB files to the destination directory.
    """
    os.makedirs(dest_dir, exist_ok=True)
    
    copied_files = []
    missing_files = []
    
    for pdb_id in pdb_ids:
        source_file = os.path.join(source_dir, f"{pdb_id}.pdb")
        dest_file = os.path.join(dest_dir, f"{pdb_id}.pdb")
        
        if os.path.exists(source_file):
            shutil.copy2(source_file, dest_file)
            copied_files.append(pdb_id)
            print(f"Copied {pdb_id}.pdb")
        else:
            missing_files.append(pdb_id)
            print(f"Missing: {pdb_id}.pdb")
    
    return copied_files, missing_files

def generate_up_files(pdb_ids, pdb_dir, up_dir):
    """
    Generate .up files from PDB files using upside_param.py
    """
    os.makedirs(up_dir, exist_ok=True)
    
    successful = []
    failed = []
    
    for pdb_id in pdb_ids:
        pdb_file = os.path.join(pdb_dir, f"{pdb_id}.pdb")
        up_file = os.path.join(up_dir, f"{pdb_id}.up")
        
        if not os.path.exists(pdb_file):
            print(f"PDB file not found: {pdb_file}")
            failed.append(pdb_id)
            continue
        
        try:
            # Run upside_param.py to generate .up file
            cmd = ["python", "upside_param.py", pdb_file, up_file]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            
            if result.returncode == 0 and os.path.exists(up_file):
                successful.append(pdb_id)
                print(f"Generated {pdb_id}.up")
            else:
                failed.append(pdb_id)
                print(f"Failed to generate {pdb_id}.up")
                if result.stderr:
                    print(f"  Error: {result.stderr.strip()}")
        
        except subprocess.TimeoutExpired:
            failed.append(pdb_id)
            print(f"Timeout generating {pdb_id}.up")
        except Exception as e:
            failed.append(pdb_id)
            print(f"Exception generating {pdb_id}.up: {e}")
    
    return successful, failed

def create_training_script(pdb_ids, training_dir):
    """
    Create a script to launch condiv training.
    """
    script_content = f"""#!/bin/bash

# Condiv training script for {len(pdb_ids)} smallest alpha proteins
# Generated automatically

echo "Starting condiv training for {len(pdb_ids)} alpha proteins"
echo "Training directory: {training_dir}"

# List of proteins
PROTEINS=({' '.join(pdb_ids)})

echo "Proteins to train:"
for protein in "${{PROTEINS[@]}}"; do
    echo "  $protein"
done

# Check if .up files exist
echo ""
echo "Checking .up files..."
missing_up=0
for protein in "${{PROTEINS[@]}}"; do
    up_file="{training_dir}/$protein.up"
    if [ -f "$up_file" ]; then
        echo "  ✓ $protein.up"
    else
        echo "  ✗ $protein.up (missing)"
        missing_up=$((missing_up + 1))
    fi
done

if [ $missing_up -gt 0 ]; then
    echo ""
    echo "Error: $missing_up .up files are missing!"
    echo "Please run the preparation script first."
    exit 1
fi

echo ""
echo "All .up files found. Ready to start training!"
echo ""

# Launch condiv training
echo "Launching condiv training..."
python condiv2.py \\
    --input-dir {training_dir} \\
    --proteins {' '.join(pdb_ids)} \\
    --output-dir results_alpha_training \\
    --epochs 1000 \\
    --batch-size 32

echo "Training completed!"
"""
    
    script_path = os.path.join(training_dir, "launch_training.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Make script executable
    os.chmod(script_path, 0o755)
    
    return script_path

def main():
    # Configuration
    selection_file = "selection/enhanced_alpha_selection_info.tsv"
    source_pdb_dir = "selection/pdbs/alpha"
    training_dir = "alpha_training_set"
    # Use all available alpha proteins
    df_all = pd.read_csv(selection_file, sep='\t')
    n_proteins = len(df_all)
    
    print(f"Preparing {n_proteins} smallest alpha proteins for condiv training")
    print("="*60)
    
    # Get smallest proteins
    print(f"Reading selection info from: {selection_file}")
    pdb_ids, df_info = get_smallest_alpha_proteins(selection_file, n_proteins)
    
    print(f"\\nSelected {len(pdb_ids)} smallest alpha proteins:")
    for i, (_, row) in enumerate(df_info.iterrows(), 1):
        pdb_id = row['PDB_Chain'].split('_')[0]
        print(f"  {i:2d}. {pdb_id} - {row['Length']} residues")
    
    # Copy PDB files
    print(f"\\nCopying PDB files to: {training_dir}")
    copied_files, missing_files = copy_pdb_files(pdb_ids, source_pdb_dir, training_dir)
    
    if missing_files:
        print(f"\\nWarning: {len(missing_files)} PDB files were missing:")
        for pdb_id in missing_files:
            print(f"  - {pdb_id}")
    
    # Generate .up files
    print(f"\\nGenerating .up files...")
    successful_up, failed_up = generate_up_files(copied_files, training_dir, training_dir)
    
    print(f"\\nSuccessfully generated {len(successful_up)} .up files")
    if failed_up:
        print(f"Failed to generate {len(failed_up)} .up files:")
        for pdb_id in failed_up:
            print(f"  - {pdb_id}")
    
    # Create training script
    print(f"\\nCreating training launch script...")
    script_path = create_training_script(successful_up, training_dir)
    print(f"Training script created: {script_path}")
    
    # Summary
    print("\\n" + "="*60)
    print("PREPARATION SUMMARY")
    print("="*60)
    print(f"Training directory: {training_dir}")
    print(f"PDB files copied: {len(copied_files)}")
    print(f".up files generated: {len(successful_up)}")
    print(f"Ready for training: {len(successful_up)} proteins")
    
    if successful_up:
        print(f"\\nTo start training, run:")
        print(f"  cd {training_dir}")
        print(f"  ./launch_training.sh")
    
    return len(successful_up) > 0

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
