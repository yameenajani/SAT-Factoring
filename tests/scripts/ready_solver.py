import sys
from ast import literal_eval
import os
import json

if __name__ == "__main__":
    bitsize = sys.argv[1]
    percentage = sys.argv[2]
    instance_number = sys.argv[3]
    include_d = int(sys.argv[4])
    
    instance_filename = "../instances/{}/original_instances/instance_{}.cnf".format(bitsize, instance_number)
    cube_filename = "../instances/{}/cubes/cubes_{}.json".format(bitsize, instance_number)
    if include_d:
        updated_instance_dir = "../instances/{}/temp_withd/{}".format(bitsize, percentage)
        updated_instance_filename = "../instances/{}/temp_withd/{}/instance_{}.cnf".format(bitsize, percentage, instance_number)
    else:
        updated_instance_dir = "../instances/{}/temp/{}".format(bitsize, percentage)
        updated_instance_filename = "../instances/{}/temp/{}/instance_{}.cnf".format(bitsize, percentage, instance_number)
    
    with open(instance_filename, "r") as ifile:
        instance = ifile.readlines()
    
    with open(cube_filename, "r") as cfile:
        data = json.load(cfile)
        cube = data['{}'.format(percentage)]

    cube = cube.split(" ")
    clauses = []
    for val in cube:
        if val != "0":
            clauses.append("\n{} 0".format(val))

    if not (os.path.exists(updated_instance_dir)):
        os.makedirs(updated_instance_dir)

    l = instance[14]
    l = l.split(" ")
    l[3] = str(int(l[3]) + len(clauses))
    instance[14] = " ".join(l) + "\n"

    with open(updated_instance_filename, "w") as file:
        for line in instance:
            file.write(line)
        for clause in clauses:
            file.write(clause)
        
