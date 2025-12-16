#!/bin/bash

# prepare_inputs.sh - Convert PDB files to condiv input format
# Usage: ./prepare_inputs.sh [pdb_folder] [output_dir]

set -e  # Exit on any error

# Default parameters
pdb_folder=${1:-"pdbs"}
output_dir=${2:-"input"}
pdb_list="pdb_list.txt"

# Check if PDB folder exists
if [ ! -d "$pdb_folder" ]; then
    echo "Error: PDB folder '$pdb_folder' not found!"
    echo "Usage: $0 [pdb_folder] [output_dir]"
    echo "Example: $0 pdbs input"
    exit 1
fi

# Setup environment (same as init.sbatch)
export MKL_THREADING_LAYER=GNU
export UPSIDE_HOME=/home/okleinmann/upside2-md/py3/
export THEANO_FLAGS='blas.ldflags='

# Activate conda environment
source /scratch/midway3/okleinmann/miniconda3/etc/profile.d/conda.sh
conda activate condiv-env

# Load modules
module load eigen/3.4
module load gcc/10.2.0
module load mpich/3.4.3+gcc-10.2.0

# Set library path
export LD_LIBRARY_PATH="/software/mpich-3.4.3-el8-x86_64+gcc-10.2.0/lib:/software/gcc-10.2.0-el8-x86_64/lib64"

# Source upside environment  
source /home/okleinmann/upside2-md/py3/source.sh

# Use conda environment Python
CONDA_PYTHON="/scratch/midway3/okleinmann/miniconda3/envs/condiv-env/bin/python"

echo "Processing PDB files from $pdb_folder..."
echo "Output directory: $output_dir"

# Create output directory
mkdir -p $output_dir

# Create pdb_list.txt header
echo "prot" > $pdb_list

# Count files
total_pdbs=$(ls -1 $pdb_folder/*.pdb 2>/dev/null | wc -l)
processed=0
successful=0

echo "Found $total_pdbs PDB files to process"
echo

# Process each PDB file
for pdb_file in $pdb_folder/*.pdb; do
    if [ -f "$pdb_file" ]; then
        # Get basename without extension
        basename_pdb=$(basename "$pdb_file" .pdb)
        processed=$((processed + 1))
        
        echo "[$processed/$total_pdbs] Converting $basename_pdb..."
        
        # Run PDB conversion
        if $CONDA_PYTHON /home/okleinmann/upside2-md/py3/py/PDB_to_initial_structure.py \
            "$pdb_file" "$output_dir/$basename_pdb" 2>/dev/null; then
            
            # Check if all required files were created
            if [ -f "$output_dir/$basename_pdb.initial.npy" ] && \
               [ -f "$output_dir/$basename_pdb.fasta" ] && \
               [ -f "$output_dir/$basename_pdb.chi" ]; then
                echo "$basename_pdb" >> $pdb_list
                successful=$((successful + 1))
                echo "  ✓ Successfully converted $basename_pdb"
            else
                echo "  ✗ Incomplete conversion for $basename_pdb (missing files)"
            fi
        else
            echo "  ✗ Failed to convert $basename_pdb (conversion error)"
        fi
    fi
done

echo
echo "=== Conversion Summary ==="
echo "Total PDB files: $total_pdbs"
echo "Processed: $processed"  
echo "Successful: $successful"
echo "Failed: $((processed - successful))"
echo
echo "Output files:"
echo "  - Input structures: $output_dir/"
echo "  - Protein list: $pdb_list"
echo
echo "Ready to run: sbatch init.sbatch"