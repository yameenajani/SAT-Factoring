#!/bin/bash

if command -v module &> /dev/null
then
    module load StdEnv/2020 gcc/9.3.0 fplll/5.4.4
    module load flint/2.9.0

    export MROOT="/home/ajaniy/scratch/coppersmith/AAAI24"
else
    export MROOT=$PWD
fi

cd core
echo "Compiling Solver"
gmake rs
echo "Compile Done"
mkdir -p ../tests/solvers/
cp -f maplesat_static ../tests/solvers/maplesat_cb
