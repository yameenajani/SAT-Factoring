#!/bin/bash

export MROOT=$PWD
cd core
echo "Compiling Solver"
gmake rs
echo "Compile Done"
cd ..
mkdir -p ../tests/solvers/
cp -f core/maplesat_static ../tests/solvers/maplesat
