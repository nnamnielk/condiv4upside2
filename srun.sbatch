#!/bin/bash
#SBATCH --job-name=condiv
#SBATCH --partition=caslake
#SBATCH --account=pi-trsosnic
#SBATCH --time=00:03:00
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=okleinmann@rcc.uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --output=condiv_training_results/slurm-%j.out
#SBATCH --error=condiv_training_results/slurm-%j.err

set -euo pipefail

# Environment
export MKL_THREADING_LAYER=GNU
export UPSIDE_HOME=/home/okleinmann/upside2-md/py3/
export THEANO_FLAGS='blas.ldflags='

# Optional: help prevent oversubscription in threaded libs
export OMP_NUM_THREADS=10

# Ensure we run from the submit directory
cd "$SLURM_SUBMIT_DIR"

# Activate conda environment with required Python packages
source ~/miniconda3/etc/profile.d/conda.sh
conda activate condiv-env

# Parameters (mirroring srun.sh)
X="$(pwd)"
Y="condiv_training_results"

mode="restart"
checkpoint="$X/$Y/initial_checkpoint.pkl"

epoch=4
batch=15              # must match #SBATCH --ntasks
n_rpx=10              # must match #SBATCH --cpus-per-task
total=$(ls selection/pdbs/alpha | wc -l)
repeat=1

# Number of restart iterations (epochs * minibatches * repeats)
step=$(( (total / batch) * epoch * repeat ))

# IMPORTANT: Do NOT wrap this with srun.
# The driver (condiv2.py) launches worker processes itself using 'srun'
# (one per protein) as job steps, consuming 1 task and n_rpx CPUs each.
python "$Y/condiv2.py" "$mode" "$checkpoint" "$step" | tee -a "$Y.output"
