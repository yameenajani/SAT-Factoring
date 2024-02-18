#!/bin/bash
while getopts b:i: flag
do
    case "${flag}" in
        b) bitsize=${OPTARG};;
        i) instance_num=${OPTARG};;
    esac
done

echo "Generating instances for ${bitsize}"
mkdir -p ../instances/$bitsize/original_instances
mkdir -p ../instances/$bitsize/cubes

python get_cnf.py ../data/$bitsize/data_$instance_num.json ../instances/$bitsize/original_instances/instance_$instance_num.cnf
python update_cnf.py ../data/$bitsize/data_$instance_num.json ../instances/$bitsize/original_instances/instance_$instance_num.cnf ../instances/$bitsize/cubes/cubes_$instance_num.json
