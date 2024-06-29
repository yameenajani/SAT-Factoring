#!/bin/bash

include_d=0
margs="-only-sat"
modedir="SAT/"
ddir=""

if [ $# -eq 0 ] || [ "$1" == "--help" ]; then
    echo "USAGE:
./run_test.sh -b {bitsize} -p {percent} -i {instance number} [-d -l]

REQUIRED ARGUMENTS:
-b {bitsize} : Specify the bitsize of the prime you want to run the test for.
-p {percent} : Specify the percentage of known bits. (Varies between 90-25% in steps of 5).
-i {instance number} : Specify the instance number you want to run.

OPTIONAL ARGUMENTS:
-d : Use the encoding with information about private exponent 'd'.
-l : Run the solver using SAT+CAS (Coppersmith called when enough low bits set)."
    exit
fi

while getopts b:p:i:dl flag
do
    case "${flag}" in
        b) bitsize=${OPTARG} ;;
        p) percent=${OPTARG} ;;
        i) instance_num=${OPTARG} ;;
        d) include_d=1; ddir="with_d/" ;;
        l) margs="-lo"; modedir="" ;;
    esac
done

if [ -z $bitsize ]; then
    echo "Need to provide bitsize of primes with -b"
    exit
fi

if [ -z $percent ]; then
    echo "Need to provide percentage of known bits with -p"
    exit
fi

if [ -z $instance_num ]; then
    echo "Need to provide instance # to solve with -i"
    exit
fi

if [ ! -f encoder/unbal/encoder ]; then
    echo "Encoder not found; compiling using GHC ..."
    cd encoder/unbal; ghc encoder; cd "$OLDPWD"
fi

data_file="data/${bitsize}/data_${instance_num}.json"
cubes_file="instances/${bitsize}/cubes/unbal/cubes_${instance_num}.json"

# Creating data file
mkdir -p data/$bitsize
cd scripts; ./gen_nums.py $bitsize $instance_num; cd "$OLDPWD"

num=$(grep n $data_file | awk '{print $2}')

# Creating cubes file
mkdir -p instances/$bitsize/cubes/unbal
./scripts/generate_cubes.py $bitsize $instance_num

printf "Generating instance ${instance_num} with ${bitsize}-bit primes"
if [ $include_d -eq 1 ]; then
    printf " (info on decryption exponent included)"
fi
printf " and setting ${percent}%% randomly chosen bits to be known\n"

encoding=$(./encoder/unbal/encoder $num n-bit recursive | ./scripts/modify_cnf.py $data_file $cubes_file $percent $include_d $instance_num)

if [ ! -f solvers/maplesat ]; then
    echo "MapleSAT not found; compiling..."
    cd ../maplesat; ./compile_maplesat.sh; cd "$OLDPWD"
fi

echo "Running MapleSAT with arguments $margs"
mkdir -p maplesat_outputs/${bitsize}/${ddir}${percent}/${modedir}
echo "$encoding" | ./solvers/maplesat $margs | tee maplesat_outputs/${bitsize}/${ddir}${percent}/${modedir}output_${instance_num}.log
