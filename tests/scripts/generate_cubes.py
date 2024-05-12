#!/usr/bin/python3

import sys
from random import sample, seed
import json
import math

def convert_to_binary(dec_val):
    return bin(dec_val).split("b")[1]

def main():
    if len(sys.argv) <= 2:
        print("Usage: {} bitsize instance_num".format(sys.argv[0]))
        quit()

    bitsize = int(sys.argv[1])
    instance_num = int(sys.argv[2])
    unbalanced = 1
    seed(instance_num)

    data_filename = "data/{}/data_{}.json".format(bitsize, instance_num)
    if unbalanced == 1:
        cubes_filename = "instances/{}/cubes/unbal/cubes_{}.json".format(bitsize, instance_num)
    else:
        cubes_filename = "instances/{}/cubes/cubes_{}.json".format(bitsize, instance_num)

    percentages = [90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25]

    cubes_obj = {}

    with open(data_filename) as fp:
        f = json.load(fp)
        p_val = f['p']
        q_val = f['q']
        n_val = f['n']

    p_list = range(bitsize,0,-1)

    if unbalanced == 1:
        output_size = math.floor(math.log(n_val, 2))+1
        q_list = range(bitsize+output_size-1,output_size-1,-1)
    else:
        q_list = range(2*bitsize,bitsize,-1)
    
    p_bin = convert_to_binary(p_val)
    q_bin = convert_to_binary(q_val)
    
    for val in percentages:
        cube = ""
        set_known_bits_count = int(len(p_list) * val/100)

        random_p_bits = sample(p_list, set_known_bits_count)
        random_q_bits = sample(q_list, set_known_bits_count)

        if p_list[-1] not in random_p_bits:
            random_p_bits.append(p_list[-1])
        if q_list[-1] not in random_q_bits:
            random_q_bits.append(q_list[-1])
        if p_list[0] not in random_p_bits:
            random_p_bits.append(p_list[0])
        if q_list[0] not in random_q_bits:
            random_q_bits.append(q_list[0])

        # Additional constraints from Heninger and Shacham
        for i in range(len(p_list)-1, -1, -1):
            if p_list[i] not in random_p_bits and q_list[i] not in random_q_bits:
                break
            elif p_list[i] in random_p_bits and q_list[i] not in random_q_bits:
                random_q_bits.append(q_list[i])
            elif p_list[i] not in random_p_bits and q_list[i] in random_q_bits:
                random_p_bits.append(p_list[i])
            else:
                continue

        for bit in random_p_bits:
            if p_bin[p_list.index(bit)] == "1":
                cube += "{} ".format(bit)
            else:
                cube += "-{} ".format(bit)
        for bit in random_q_bits:
            if q_bin[q_list.index(bit)] == "1":
                cube += "{} ".format(bit)
            else:
                cube += "-{} ".format(bit)

        # Set high bits of p to 0 if unbalanced encoding is used
        if unbalanced == 1:
            for bit in range(bitsize+1,output_size):
                cube += "-{} ".format(bit)

        cube += "0"

        cubes_obj.update({"{}".format(val): cube})

    with open(cubes_filename, "w") as cube_f:
        cube_f.write(json.dumps(cubes_obj, indent=4))

if __name__ == '__main__':
    main()
