#!/bin/bash
while getopts b:n: flag
do
    case "${flag}" in
        b) bitsize=${OPTARG};;
        n) num_instances=${OPTARG};;
    esac
done

echo "Generating instances for ${bitsize}"
mkdir -p ../instances/$bitsize/new_instances
mkdir -p ../instances/$bitsize/cubes
for ((i=1; i<=$num_instances; i++))
do
    python get_cnf.py ../data/$bitsize/data_$i.json ../instances/$bitsize/new_instances/instance_$i.cnf
    python update_cnf.py ../data/$bitsize/data_$i.json ../instances/$bitsize/new_instances/instance_$i.cnf ../instances/$bitsize/cubes/cubes_$i.json
done
    # done
# done
