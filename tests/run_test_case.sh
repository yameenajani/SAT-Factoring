#!/bin/bash
include_d=0
while getopts b:p:i:d flag
do
    case "${flag}" in
        b) bitsize=${OPTARG};;
        p) percent=${OPTARG};;
        i) instance=${OPTARG};;
        d) include_d=1;;
    esac
done

if [ -z $bitsize ]
then
    echo "Need to provide bitsize of primes with -b"
    exit
fi

if [ -z $percent ]
then
    echo "Need to provide percentage of known bits with -p"
    exit
fi

if [ -z $instance ]
then
    echo "Need to provide instance # to solve with -i"
    exit
fi

cd scripts
if [ $include_d -eq 1 ]
then
    command="./init_solver.sh -b $bitsize -p $percent -i $instance -d"
else
    command="./init_solver.sh -b $bitsize -p $percent -i $instance"
fi
echo $command
eval $command
cd ..

mkdir -p outputs
mkdir -p outputs/$bitsize
if [ $include_d -eq 1 ]
then
    mkdir -p outputs/$bitsize/with_d/$percent
    mkdir -p outputs/$bitsize/with_d/$percent/SAT
else
    mkdir -p outputs/$bitsize/$percent
    mkdir -p outputs/$bitsize/$percent/SAT
fi

echo "Running only SAT"
if [ $include_d -eq 1 ]
then
    command="./solvers/maplesat_cb -only-sat instances/$bitsize/temp_withd/$percent/instance_$instance.cnf | tee outputs/$bitsize/with_d/$percent/SAT/output_$instance.log"
else
    command="./solvers/maplesat_cb -only-sat instances/$bitsize/temp/$percent/instance_$instance.cnf | tee outputs/$bitsize/$percent/SAT/output_$instance.log"
fi
echo $command
eval $command
echo "Finished running only SAT"

echo "Running SAT+CAS"
if [ $include_d -eq 1 ]
then
    command="./solvers/maplesat_cb instances/$bitsize/temp_withd/$percent/instance_$instance.cnf | tee outputs/$bitsize/with_d/$percent/output_$instance.log"
else
    command="./solvers/maplesat_cb instances/$bitsize/temp/$percent/instance_$instance.cnf | tee outputs/$bitsize/$percent/output_$instance.log"
fi
echo $command
eval $command
echo "Finished running SAT+CAS"

echo "Output files available in tests/outputs/$bitsize/ dir"
