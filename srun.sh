#!/bin/bash
#SBATCH --job-name=condiv_run
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=6
#SBATCH --time=35:00:00
#SBATCH --mail-user=okleinmann@rcc.uchicago.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=condiv_run_%j.out
#SBATCH --error=condiv_run_%j.err

export MKL_THREADING_LAYER=GNU

X=`pwd`
Y=test_00

mode=restart
checkpoint=$X/$Y/initial_checkpoint.pkl

batch_size=6 # number of proteins total
n_rpx=6 # number of cpu cores
step=5 # number of training loops

echo python3.7 $X/$Y/condiv2.py $mode $checkpoint $step | tee -a $X/$Y.output
python3.7 $X/$Y/condiv2.py $mode $checkpoint $step | tee -a $X/$Y.output
