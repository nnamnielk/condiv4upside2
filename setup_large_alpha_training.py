#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import shutil
from Bio.PDB import PDBParser

def get_alpha_proteins(selection_file, n=50):
    """
    Get the n smallest alpha proteins from the selection info file.
    """
    df = pd.read_csv(selection_file, sep='\t')
    
    # Sort by length and get the smallest n proteins
    df_sorted = df.sort_values('Length').head(n)
    
    # Extract unique PDB IDs
    pdb_ids = []
    seen = set()
    for pdb_chain in df_sorted['PDB_Chain']:
        pdb_id = pdb_chain.split('_')[0]
        if pdb_id not in seen:
            pdb_ids.append(pdb_id)
            seen.add(pdb_id)
    
    return pdb_ids, df_sorted

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

def setup_training_proteins(pdb_ids, source_pdb_dir, input_dir):
    """
    Set up all training proteins with required files
    """
    os.makedirs(input_dir, exist_ok=True)
    
    successful_proteins = []
    
    for pdb_id in pdb_ids:
        try:
            pdb_file = os.path.join(source_pdb_dir, f"{pdb_id}.pdb")
            if not os.path.exists(pdb_file):
                print(f"Warning: {pdb_file} not found, skipping")
                continue
            
            # Create .initial.npy file in input/ directory
            initial_file = os.path.join(input_dir, f"{pdb_id}.initial.npy")
            n_res = pdb_to_initial_npy(pdb_file, initial_file)
            
            successful_proteins.append(pdb_id)
            print(f"✓ Prepared {pdb_id}: {n_res} residues")
            
        except Exception as e:
            print(f"✗ Failed to prepare {pdb_id}: {e}")
            continue
    
    return successful_proteins

def update_condiv_parameters(minibatch_size=15, n_threads=10):
    """
    Update condiv2.py with new parameters
    """
    with open('condiv2.py', 'r') as f:
        content = f.read()
    
    # Update minibatch_size
    content = content.replace(
        'minibatch_size = 11',
        f'minibatch_size = {minibatch_size}'
    )
    
    # Update n_threads
    content = content.replace(
        'n_threads = 6 # this is per protein, must be >3',
        f'n_threads = {n_threads} # this is per protein, must be >3'
    )
    
    with open('condiv2.py', 'w') as f:
        f.write(content)
    
    print(f"Updated condiv2.py: minibatch_size={minibatch_size}, n_threads={n_threads}")

def update_slurm_resources(minibatch_size, n_threads):
    """
    Update srun.sh with correct CPU allocation
    """
    total_cpus = minibatch_size * n_threads
    
    with open('srun.sh', 'r') as f:
        content = f.read()
    
    # Update ntasks
    content = content.replace(
        '--ntasks=66',
        f'--ntasks={total_cpus}'
    )
    
    with open('srun.sh', 'w') as f:
        f.write(content)
    
    print(f"Updated srun.sh: --ntasks={total_cpus} (for {minibatch_size} proteins × {n_threads} threads)")

def create_pdb_list(protein_ids):
    """
    Create pdb_list file
    """
    with open('pdb_list', 'w') as f:
        for pdb_id in protein_ids:
            f.write(f"{pdb_id}\n")
    
    print(f"Created pdb_list with {len(protein_ids)} proteins")

def main():
    # Configuration
    n_proteins = 50  # Number of proteins to use
    minibatch_size = 15
    n_threads = 10
    
    print(f"Setting up large alpha training set")
    print(f"Proteins: {n_proteins}")
    print(f"Minibatch size: {minibatch_size}")
    print(f"Threads per protein: {n_threads}")
    print("=" * 60)
    
    # Get alpha proteins
    selection_file = "selection/enhanced_alpha_selection_info.tsv"
    source_pdb_dir = "selection/pdbs/alpha"
    input_dir = "input"
    
    print(f"Reading selection info from: {selection_file}")
    pdb_ids, df_info = get_alpha_proteins(selection_file, n_proteins)
    
    print(f"\nSelected {len(pdb_ids)} unique alpha proteins (smallest by length)")
    
    # Setup training proteins
    print(f"\nPreparing proteins in: {input_dir}")
    successful_proteins = setup_training_proteins(pdb_ids, source_pdb_dir, input_dir)
    
    if not successful_proteins:
        print("Error: No proteins were successfully prepared!")
        return False
    
    print(f"\nSuccessfully prepared {len(successful_proteins)} proteins")
    
    # Create pdb_list
    create_pdb_list(successful_proteins)
    
    # Update condiv2.py parameters
    update_condiv_parameters(minibatch_size, n_threads)
    
    # Update SLURM resources
    update_slurm_resources(minibatch_size, n_threads)
    
    # Summary
    print("\n" + "=" * 60)
    print("SETUP COMPLETE")
    print("=" * 60)
    print(f"Proteins prepared: {len(successful_proteins)}")
    print(f"Minibatch size: {minibatch_size}")
    print(f"Threads per protein: {n_threads}")
    print(f"Total CPU cores needed: {minibatch_size * n_threads}")
    print(f"Input directory: {input_dir}")
    print(f"PDB list: pdb_list")
    
    print("\nTo start training:")
    print("  sbatch init.sh    # Initialize training")
    print("  sbatch srun.sh    # Start training (after init completes)")
    
    return True

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
