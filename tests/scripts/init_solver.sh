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

if [ $include_d -eq 1 ]
then
    echo "Including d"
    python ready_solver.py $bitsize $percent $instance 1
    python add_constraints.py $bitsize $percent $instance
else
    python ready_solver.py $bitsize $percent $instance 0
fi


