# Upside 2 ConDiv Training

This project trains the Upside Force Field using the dual-target contrastive divergence (ConDiv) method described in Dr. Peng's 2021 ACS paper.

***

## Source Code
The full source code for this project is available on GitHub.
- **GitHub Repository:** [nnamnielk/condiv4upside2](https://github.com/nnamnielk/condiv4upside2/tree/main)

***

## Prerequisites

### Project Structure
The project directory structure:

```
condiv/
├── init_param/               # Initial force field parameters
│   ├── environment.h5
│   ├── bb_env.dat
│   ├── sidechain.h5
│   ├── hbond.h5
│   └── sheet
├── input/                    # Training input files (generated from PDBs)
│   ├── <pdbid>.initial.npy  # Initial structure coordinates
│   ├── <pdbid>.fasta        # Sequence files
│   └── <pdbid>.chi          # Chi angles (side chain rotamers)
├── training/                 # Training set preparation (optional)
│   └── pengxd/400_700/      # Example: all-alpha proteins 400-700 residues
│       ├── pdb/             # Source PDB files
│       ├── input/           # Generated input files
│       └── pdb_list.txt     # List of proteins in this set
├── pdb_list.txt             # Active protein list for training
├── condiv2.py               # Main training script
├── prepare_inputs.sh        # PDB to input conversion script
├── init.sbatch              # Initialize training checkpoint
└── srun.sbatch              # Run training iterations
```

### Environment
- The `$UPSIDE_HOME` environment variable must be set to your Upside installation
- Location: `/home/okleinmann/upside2-md/py3/`
- You may need to update the `py/` directory with the contents of `updated_utils/`

***

## Step 1: Installation

1.  Create the Conda environment:
    ```bash
    conda env create -f environment/conda-env2.yml
    ```
2.  Activate the new environment:
    ```bash
    conda activate <your-env-name>
    ```
3.  **Find and replace** all hardcoded references to `/home/okleinmann` in the project files.
4.  If not using a Slurm workload manager, remove all `sbatch` and `srun` commands from the scripts.

***

## Step 2: Training Workflow

Follow these steps to prepare the data and run the training.

### 1. Select a Training Set
To generate a training set using SCOPe, run the Jupyter Notebook from within your Conda environment:
```bash
jupyter notebook notebooks/1.training-selection.ipynb
```
This will help you select proteins and download PDB files.

### 2. Prepare Input Files
Convert PDB files to Upside input format using `prepare_inputs.sh`:

```bash
./prepare_inputs.sh [pdb_folder] [output_dir]
# Example: ./prepare_inputs.sh training/pengxd/400_700/pdb input
```

```bash
./pdb_list_gen.sh [pdb_folder] > pdb_list.txt
```

This script:
- Uses `/home/okleinmann/upside2-md/py3/py/PDB_to_initial_structure.py` to convert each PDB
- Creates three files per protein: `.initial.npy`, `.fasta`, `.chi`
- Generates `pdb_list.txt` with successfully converted proteins
- Outputs to `input/` directory (or specified output directory)

**Note:** The training scripts (`init.sbatch` and `srun.sbatch`) reference files in the working directory:
- `upside_input_dir=input` (line 23 in init.sbatch)
- `pdb_list=pdb_list.txt` (line 24 in init.sbatch)

If you prepare inputs in a subdirectory like `training/pengxd/400_700/`, copy the files to the working directory:
```bash
cp training/pengxd/400_700/pdb_list.txt .
cp -r training/pengxd/400_700/input .
```

### 3. Initialize Training Checkpoint
Generate the initial checkpoint file for training:
```bash
sbatch init.sbatch
```

This script:
- Reads proteins from `input/` directory and `pdb_list.txt`
- Loads initial parameters from `init_param/`
- Creates `condiv_training_results/initial_checkpoint.pkl`

Mode: `initialize`

### 4. Run the Training
Execute the dual-target ConDiv training:
```bash
sbatch srun.sbatch
```

This script:
- Loads the checkpoint from step 3
- Runs training iterations (epochs × minibatches × repeats)
- Launches worker processes using `srun` for parallel protein simulations
- Updates force field parameters after each minibatch

Mode: `restart`

***

## Key File Locations

### Training Directory Reference
- `training/pengxd/400_700/` - Example training set (not directly used by code)
- Training scripts use local copies: `input/` and `pdb_list.txt` in working directory

### PDB Conversion
- Script: `$UPSIDE_HOME/py3/py/PDB_to_initial_structure.py`#you must be running Kleinmann's fork for the correct version
- Called by: `prepare_inputs.sh`
- Generates: `.initial.npy`, `.fasta`, `.chi` files

### Configuration
`init.sbatch` specifies:
- `upside_input_dir=input` (line 23)
- `pdb_list=pdb_list.txt` (line 24)
- `initial_param=init_param` (line 22)

***

## Step 3: Analyze Results
After training is complete, analyze the loss to evaluate the model's performance.

***

## Citation
Peng, X., et al. (2021). Prediction and Validation of a Protein’s Free Energy Surface Using Hydrogen Exchange and (Importantly) Its Denaturant Dependence. *Journal of Chemical Theory and Computation*, 17(11), 7052–7062.
[https://doi.org/10.1021/acs.jctc.1c00960](https://doi.org/10.1021/acs.jctc.1c00960)
