#!/bin/bash

./gcm.e > gcm.out 2>gcm.err
cp restart.nc start.nc
cp restartfi.nc startfi.nc
nfiles=$(ls -l snapshot*.npy | wc -l)
python snapshot.py diagfi.nc
mv gcm.out gcm$nfiles.out
mv gcm.err gcm$nfiles.err
mv snapshot.npy snapshot$nfiles.npy
cp snapshot$nfiles.npy ../output/$1_snapshot.npy
if [ $nfiles -lt $2 ]
then
        cp -a * $PBS_O_WORKDIR/
        rm -rf *
        cd $PBS_O_WORKDIR/
	qsub runlmdz
else
        cp -a * $PBS_O_WORKDIR/
        rm -rf *
        cd $PBS_O_WORKDIR/
	python release.py
fi
