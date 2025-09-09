#!/bin/bash
#SBATCH --job-name=condiv
#SBATCH --partition=caslake
#SBATCH --account=pi-trsosnic
#SBATCH --mem-per-cpu=3G
#SBATCH --time=01:00:00
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=okleinmann@rcc.uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --output=condiv_training_results/slurm-%j.out
#SBATCH --error=condiv_training_results/slurm-%j.err

set -eo pipefail 

# Environment
export MKL_THREADING_LAYER=GNU
export UPSIDE_HOME=/home/okleinmann/upside2-md/py3/
export THEANO_FLAGS='blas.ldflags='

# Expand SLURM port range to use full available range
export SLURM_SRUN_COMM_PORT_RANGE=0-65535

# Optional: help prevent oversubscription in threaded libs
export OMP_NUM_THREADS=10

# Ensure we run from the submit directory
cd "$SLURM_SUBMIT_DIR"

# Activate conda environment first
source /scratch/midway3/okleinmann/miniconda3/etc/profile.d/conda.sh
conda activate condiv-env

# Then load modules (this will override conda's library paths)
module load eigen/3.4
module load gcc/10.2.0
module load mpich/3.4.3+gcc-10.2.0

# Set library path to prioritize the correct MPI libraries AFTER conda activation
# Exclude conda lib to avoid GLIBCXX conflicts
export LD_LIBRARY_PATH="/software/mpich-3.4.3-el8-x86_64+gcc-10.2.0/lib:/software/gcc-10.2.0-el8-x86_64/lib64"

# Source upside environment  
source /home/okleinmann/upside2-md/py3/source.sh

# Ensure srun exports all environment variables to worker nodes
export SLURM_EXPORT_ENV=ALL



# Parameters (mirroring srun.sh)
X="$(pwd)"
Y="condiv_training_results"

mode="restart"
checkpoint="$X/$Y/initial_checkpoint.pkl"

epoch=4
batch=2              # must match #SBATCH --ntasks
n_rpx=4              # must match #SBATCH --cpus-per-task
total=405
repeat=1

# Number of restart iterations (epochs * minibatches * repeats)
step=$(( (total / batch) * epoch * repeat ))

# Use conda environment Python directly
CONDA_PYTHON="/scratch/midway3/okleinmann/miniconda3/envs/condiv-env/bin/python"

# IMPORTANT: Do NOT wrap this with srun.
# The driver (condiv2.py) launches worker processes itself using 'srun'
# (one per protein) as job steps, consuming 1 task and n_rpx CPUs each.
$CONDA_PYTHON "$Y/condiv2.py" "$mode" "$checkpoint" "$step" | tee -a "$Y.output"
