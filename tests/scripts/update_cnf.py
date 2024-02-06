from ast import literal_eval
import linecache
import sys
from random import randint, sample
import json

def convert_to_binary(dec_val):
    return bin(dec_val).split("b")[1]

def main():
    input_filename = sys.argv[1]
    second_filename = sys.argv[2]
    cubes_filename = sys.argv[3]

    percentages = [90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25]

    clauses = []
    cubes_obj = {}

    q_list = linecache.getline(second_filename, 7)
    q_list = q_list.split(": ")
    q_list = literal_eval(q_list[1])

    p_list = linecache.getline(second_filename, 6)
    p_list = p_list.split(": ")
    p_list = literal_eval(p_list[1])

    make_msb_0_len = len(p_list) - len(q_list)

    for i in range(0, make_msb_0_len):
        clauses.append("\n-{} 0".format(p_list[i]))

    p_list = p_list[make_msb_0_len:]

    with open(input_filename) as fp:
        f = json.load(fp)
        p_val = f['p']
        q_val = f['q']
        n_val = f['n']

    comments = []
    comments.append("c {}\n".format(n_val))
    comments.append("c {}\n".format(p_list[-1]))
    comments.append("c {}\n".format(p_list[0]))
    comments.append("c {}\n".format(q_list[-1]))
    comments.append("c {}\n".format(q_list[0]))
    
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

        cube += "0"

        cubes_obj.update({"{}".format(val): cube})

    with open(cubes_filename, "w") as cube_f:
        cube_f.write(json.dumps(cubes_obj, indent=4))

    with open(second_filename, "a") as f:
        for clause in clauses:
            f.write(clause)

    with open(second_filename, "r") as file:
        lines = file.readlines()

    l = lines[9]
    l = l.split(" ")
    l[3] = str(int(l[3]) + len(clauses))
    lines[9] = " ".join(l) + "\n"

    with open(second_filename, "w") as file:
        for c in comments:
            file.write(c)
        for line in lines:
            file.write(line)

if __name__ == '__main__':
    main()