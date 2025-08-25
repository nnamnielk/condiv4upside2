#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import shutil
from Bio.PDB import PDBParser
import subprocess

def pdb_to_initial_npy(pdb_file, output_file):
    """
    Convert PDB file to .initial.npy format expected by condiv2.py
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    # Extract CA coordinates
    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id('CA'):
                    ca_atom = residue['CA']
                    ca_coords.append(ca_atom.get_coord())
    
    if not ca_coords:
        raise ValueError(f"No CA atoms found in {pdb_file}")
    
    # Convert to numpy array and save
    coords_array = np.array(ca_coords, dtype=np.float32)
    np.save(output_file, coords_array)
    
    return len(ca_coords)

def create_fasta_from_pdb(pdb_file, fasta_file):
    """
    Create FASTA file from PDB structure
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    # Standard amino acid codes
    aa_codes = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    sequence = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == ' ':  # standard residue
                    res_name = residue.get_resname()
                    if res_name in aa_codes:
                        sequence.append(aa_codes[res_name])
                    else:
                        sequence.append('X')  # unknown residue
    
    if not sequence:
        raise ValueError(f"No valid residues found in {pdb_file}")
    
    # Write FASTA file
    with open(fasta_file, 'w') as f:
        f.write(f">{os.path.basename(pdb_file).replace('.pdb', '')}\n")
        f.write(''.join(sequence) + '\n')
    
    return ''.join(sequence)

def create_chi_file(pdb_id, chi_file):
    """
    Create a dummy chi file (side chain angles)
    """
    with open(chi_file, 'w') as f:
        f.write(f"# Chi angles for {pdb_id}\n")
        f.write("# Dummy file - no chi angles specified\n")

def setup_training_proteins(training_dir, pdb_ids):
    """
    Set up all training proteins with required files
    """
    os.makedirs(training_dir, exist_ok=True)
    
    successful_proteins = []
    protein_info = []
    
    for pdb_id in pdb_ids:
        try:
            pdb_file = f"alpha_training_set/{pdb_id}.pdb"
            if not os.path.exists(pdb_file):
                print(f"Warning: {pdb_file} not found, skipping")
                continue
            
            # Create required files
            base_name = os.path.join(training_dir, pdb_id)
            initial_file = f"{base_name}.initial.npy"
            fasta_file = f"{base_name}.fasta"
            chi_file = f"{base_name}.chi"
            
            # Convert PDB to initial coordinates
            n_res = pdb_to_initial_npy(pdb_file, initial_file)
            
            # Create FASTA sequence
            sequence = create_fasta_from_pdb(pdb_file, fasta_file)
            
            # Create chi file
            create_chi_file(pdb_id, chi_file)
            
            successful_proteins.append(pdb_id)
            protein_info.append({
                'pdb_id': pdb_id,
                'n_residues': n_res,
                'sequence_length': len(sequence),
                'sequence': sequence
            })
            
            print(f"✓ Prepared {pdb_id}: {n_res} residues")
            
        except Exception as e:
            print(f"✗ Failed to prepare {pdb_id}: {e}")
            continue
    
    return successful_proteins, protein_info

def create_protein_list(training_dir, protein_ids):
    """
    Create protein list file for condiv2.py
    """
    list_file = os.path.join(training_dir, 'protein_list.txt')
    with open(list_file, 'w') as f:
        f.write("prot\n")  # header required by condiv2.py
        for pdb_id in protein_ids:
            f.write(f"{pdb_id}\n")
    return list_file

def create_launch_script(training_dir, protein_list_file, init_param_dir):
    """
    Create script to launch condiv training
    """
    script_content = f"""#!/bin/bash

# Condiv training launch script
# Generated automatically

echo "Starting condiv training setup..."

# Set environment variables
export UPSIDE_HOME="/project2/trsosnic/upside2"
export PYTHONPATH="$UPSIDE_HOME/py:$PYTHONPATH"

# Check if UPSIDE_HOME exists
if [ ! -d "$UPSIDE_HOME" ]; then
    echo "Error: UPSIDE_HOME directory not found: $UPSIDE_HOME"
    echo "Please set the correct path to your Upside installation"
    exit 1
fi

# Parameters
INIT_DIR="{init_param_dir}"
PROTEIN_DIR="{training_dir}"
PROTEIN_LIST="{protein_list_file}"
OUTPUT_DIR="condiv_training_results"

echo "Configuration:"
echo "  Init parameters: $INIT_DIR"
echo "  Protein directory: $PROTEIN_DIR"
echo "  Protein list: $PROTEIN_LIST"
echo "  Output directory: $OUTPUT_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Initialize condiv training
echo "Initializing condiv training..."
python condiv2.py initialize "$INIT_DIR" "$PROTEIN_DIR" "$PROTEIN_LIST" "$OUTPUT_DIR"

if [ $? -eq 0 ]; then
    echo "Initialization successful!"
    echo ""
    echo "To start training, run:"
    echo "  python condiv2.py restart $OUTPUT_DIR/initial_checkpoint.pkl 100"
    echo ""
    echo "Training files are in: $OUTPUT_DIR"
else
    echo "Initialization failed!"
    exit 1
fi
"""
    
    script_path = os.path.join(training_dir, "launch_condiv.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    os.chmod(script_path, 0o755)
    return script_path

def main():
    print("Setting up condiv training for alpha proteins")
    print("=" * 50)
    
    # Get the smallest alpha proteins
    selection_file = "selection/enhanced_alpha_selection_info.tsv"
    df = pd.read_csv(selection_file, sep='\t')
    # Use all alpha proteins from selection
    df_sorted = df.sort_values('Length')
    
    pdb_ids = []
    for pdb_chain in df_sorted['PDB_Chain']:
        pdb_id = pdb_chain.split('_')[0]
        pdb_ids.append(pdb_id)
    
    # Remove duplicates while preserving order
    unique_pdb_ids = []
    seen = set()
    for pdb_id in pdb_ids:
        if pdb_id not in seen:
            unique_pdb_ids.append(pdb_id)
            seen.add(pdb_id)
    
    print(f"Selected {len(unique_pdb_ids)} unique proteins:")
    for i, pdb_id in enumerate(unique_pdb_ids, 1):
        length = df[df['PDB_Chain'].str.startswith(pdb_id)]['Length'].iloc[0]
        print(f"  {i:2d}. {pdb_id} ({length} residues)")
    
    # Setup training directory
    training_dir = "condiv_training_proteins"
    successful_proteins, protein_info = setup_training_proteins(training_dir, unique_pdb_ids)
    
    if not successful_proteins:
        print("Error: No proteins were successfully prepared!")
        return False
    
    print(f"\nSuccessfully prepared {len(successful_proteins)} proteins")
    
    # Create protein list file
    protein_list_file = create_protein_list(training_dir, successful_proteins)
    print(f"Created protein list: {protein_list_file}")
    
    # Create launch script
    init_param_dir = "init_param"
    launch_script = create_launch_script(training_dir, protein_list_file, init_param_dir)
    print(f"Created launch script: {launch_script}")
    
    # Summary
    print("\n" + "=" * 50)
    print("SETUP COMPLETE")
    print("=" * 50)
    print(f"Training directory: {training_dir}")
    print(f"Proteins prepared: {len(successful_proteins)}")
    print(f"Launch script: {launch_script}")
    
    print("\nTo start training:")
    print(f"  cd {training_dir}")
    print(f"  ./launch_condiv.sh")
    
    return True

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
