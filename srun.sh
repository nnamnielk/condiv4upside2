export MKL_THREADING_LAYER=GNU
export UPSIDE_HOME=/home/okleinmann/upside2-md/py3/
export THEANO_FLAGS='blas.ldflags='

X=pwd
Y=condiv_training_results

mode=restart
checkpoint=$X/$Y/initial_checkpoint.pkl

epoch=4
batch=15
n_rpx=10
total=$(ls selection/pdbs/alpha | wc -l)
repeat=1

p=caslake

let step=(total/batch)*epoch*repeat

#caslake
salloc --partition $p -t 00:10:00 --account=pi-trsosnic --job-name=condiv --ntasks=$batch --cpus-per-task=$n_rpx python $Y/condiv2.py $mode $checkpoint $step | tee -a $Y.output


