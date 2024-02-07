#!/bin/bash
include_d=0
while getopts b:p:d: flag
do
    case "${flag}" in
        b) bitsize=${OPTARG};;
        p) percent=${OPTARG};;
        d) include_d=1;;
    esac
done

cd scripts
python gen_nums.py $bitsize 1
./get_instance.sh -b $bitsize -n 1
if [ $include_d -eq 1 ]
then
    ./init_solver.sh -b $bitsize -p $percent -i 1 -d true
else
    ./init_solver.sh -b $bitsize -p $percent -i 1
fi
cd ../
echo "Running SAT+CAS"
if [ $include_d -eq 1 ]
then
    ./solvers/maplesat_cb instances/$bitsize/temp_withd/$percent/instance_1.cnf > outputs/$bitsize/with_d/$percent/output_1.log
else
    ./solvers/maplesat_cb instances/$bitsize/temp/$percent/instance_1.cnf > outputs/$bitsize/$percent/output_1.log
fi
echo "Finished running SAT+CAS"
echo "Running only SAT"
if [ $include_d -eq 1 ]
then
    ./solvers/maplesat_cb -only-sat instances/$bitsize/temp_withd/$percent/instance_1.cnf > outputs/$bitsize/with_d/$percent/SAT/output_1.log
else
    ./solvers/maplesat_cb -only-sat instances/$bitsize/temp/$percent/instance_1.cnf > outputs/$bitsize/$percent/SAT/output_1.log
fi
echo "Finished running only SAT"
echo "Output files available in tests/outputs/$bitsize/ dir"
