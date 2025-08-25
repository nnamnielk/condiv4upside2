#!/bin/bash
export MKL_THREADING_LAYER=GNU
export UPSIDE_HOME=/home/okleinmann/upside2-md/py3/
export THEANO_FLAGS='blas.ldflags='
python /home/okleinmann/upside2-ipy/condiv/condiv_training_results/ConDiv.py restart /home/okleinmann/upside2-ipy/condiv/condiv_training_results/initial_checkpoint.pkl 36
