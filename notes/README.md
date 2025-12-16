# Upside 2 ConDiv Training

This project trains the Upside Force Field using the dual-target contrastive divergence (ConDiv) method described in Dr. Peng's 2021 ACS paper.

***

## Source Code
The full source code for this project is available on GitHub.
- **GitHub Repository:** [nnamnielk/condiv4upside2](https://github.com/nnamnielk/condiv4upside2/tree/main)

***

## Prerequisites

### Project Structure
A minimal project directory will look like this:

```
upside-2-condiv/
├── init_param/
│   └── (Parameters for initialization)
├── training/
│   └── PDB/
│       ├── protein1.pdb
│       ├── protein2.pdb
│       └── ...
├── condiv2.py
├── init.sbatch
└── srun.sbatch
```
### Environment
- The `$UPSIDE_HOME` environment variable must be set.
- You may need to update the `py/` directory with the contents of `updated_utils/`.

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
- `notebook/1.training-selection.ipynb`

### 2. Prepare Inputs
This script creates `.chi`, `.fasta`, and `.initial.npy` files in the `input/` directory. It also generates a PDB manifest (`pdb_list.txt`).
- `prepare_inputs.sh`

### 3. Initialize Upside Files
This script generates the `.up` files required for each minibatch.
- `init.sbatch`

### 4. Run the Training
This script executes the dual-target ConDiv training.
- `srun.sbatch`

***

## Step 3: Analyze Results
After training is complete, analyze the loss to evaluate the model's performance.

***

## Citation
Peng, X., et al. (2021). Prediction and Validation of a Protein’s Free Energy Surface Using Hydrogen Exchange and (Importantly) Its Denaturant Dependence. *Journal of Chemical Theory and Computation*, 17(11), 7052–7062.
[https://doi.org/10.1021/acs.jctc.1c00960](https://doi.org/10.1021/acs.jctc.1c00960)
