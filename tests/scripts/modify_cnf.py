#!/usr/bin/python3

import sys
import json
from ast import literal_eval
from random import sample, seed

# fetch_data reads the encoding and fetches the lists of variables corresponding to p, q and d
def fetch_data(encoding):
    q_list = encoding[11]
    q_list = q_list.split(": ")
    q_list = literal_eval(q_list[1])

    p_list = encoding[10]
    p_list = p_list.split(": ")
    p_list = literal_eval(p_list[1])

    N_size = encoding[8]
    N_size = N_size.split(" ")
    N = int(N_size[5])
    N_size = literal_eval(N_size[6])
    N_size = len(N_size)

    num_variables = encoding[14]
    num_variables = num_variables.split(" ")
    num_variables = int(num_variables[2])

    d_list = [*range(num_variables+1, num_variables+N_size+1)]
    d_list.reverse()

    num_variables += N_size

    return p_list, q_list, d_list, num_variables, N

# and2 adds clauses corresponding to the logical expression C = A and B
def and2(clauses, out, in_1, in_2):
    clauses.append("\n{} -{} -{} 0".format(out, in_1, in_2))
    clauses.append("\n-{} {} 0".format(out, in_1))
    clauses.append("\n-{} {} 0".format(out, in_2))

# or2 adds clauses corresponding to the logical expression C = A or B
def or2(clauses, out, in_1, in_2):
    clauses.append("\n-{} {} {} 0".format(out, in_1, in_2))
    clauses.append("\n{} -{} 0".format(out, in_1))
    clauses.append("\n{} -{} 0".format(out, in_2))

# xor2 adds clauses corresponding to the logical expression C = A xor B
def xor2(clauses, out, in_1, in_2):
    clauses.append("\n-{} -{} -{} 0".format(out, in_1, in_2))
    clauses.append("\n{} -{} {} 0".format(out, in_1, in_2))
    clauses.append("\n{} {} -{} 0".format(out, in_1, in_2))
    clauses.append("\n-{} {} {} 0".format(out, in_1, in_2))

# xor3 adds clauses corresponding to the logical expression D = A xor B xor C
def xor3(clauses, out, in_1, in_2, in_3):
    clauses.append("\n{} -{} -{} -{} 0".format(out, in_1, in_2, in_3))
    clauses.append("\n-{} -{} -{} {} 0".format(out, in_1, in_2, in_3))
    clauses.append("\n-{} -{} {} -{} 0".format(out, in_1, in_2, in_3))
    clauses.append("\n{} -{} {} {} 0".format(out, in_1, in_2, in_3))
    clauses.append("\n-{} {} -{} -{} 0".format(out, in_1, in_2, in_3))
    clauses.append("\n{} {} -{} {} 0".format(out, in_1, in_2, in_3))
    clauses.append("\n{} {} {} -{} 0".format(out, in_1, in_2, in_3))
    clauses.append("\n-{} {} {} {} 0".format(out, in_1, in_2, in_3))

# exp3 adds clauses corresponding to the logical expression D = (A or B) and (A or C) and (B or C)
def exp3(clauses, out, in_1, in_2, in_3):
    clauses.append("\n-{} {} {} 0".format(out, in_1, in_2))
    clauses.append("\n-{} {} {} 0".format(out, in_1, in_3))
    clauses.append("\n-{} {} {} 0".format(out, in_2, in_3))
    clauses.append("\n{} -{} -{} 0".format(out, in_2, in_3))
    clauses.append("\n{} -{} -{} 0".format(out, in_1, in_3))
    clauses.append("\n{} -{} -{} 0".format(out, in_1, in_2))

# Half Adder (Adds 2 bits)
def half_adder(clauses, s, c, in_1, in_2):
    xor2(clauses, s, in_1, in_2)
    and2(clauses, c, in_1, in_2)

# Full Adder (Adds 3 bits)
def full_adder(clauses, s, c, in_1, in_2, in_3):
    xor3(clauses, s, in_1, in_2, in_3)
    exp3(clauses, c, in_1, in_2, in_3)

def add_vars(clauses, vars, add_list):
    for i in range(len(add_list)):
        l = add_list[i]
        
        if len(l) >= 3:
            x1 = l.pop(0)
            x2 = l.pop(0)
            x3 = l.pop(0)

            s = vars+1
            c = vars+2
            
            full_adder(clauses, s, c, x1, x2, x3)

            add_list[i].append(s)
            if i+1 >= len(add_list):
                add_list.append([c])
            else:
                add_list[i+1].append(c)

            vars += 2

            return clauses, vars, False
        
        if len(l) >= 2:
            x1 = l.pop(0)
            x2 = l.pop(0)

            s = vars + 1
            c = vars + 2

            half_adder(clauses, s, c, x1, x2)

            add_list[i].append(s)
            if i+1 >= len(add_list):
                add_list.append([c])
            else:
                add_list[i+1].append(c)

            vars += 2

            return clauses, vars, False

    return clauses, vars, True

def convert_to_binary(dec_val, width):
    return bin(dec_val).split("b")[1].zfill(width)

if __name__ == "__main__":
    if len(sys.argv) <= 5:
        print("Usage: {} data_filename cube_filename percentage include_d instance_num".format(sys.argv[0]))
        print("The original unmodified instance is provided via the standard input")
        quit()

    data_filename = sys.argv[1]
    cube_filename = sys.argv[2]
    percentage = sys.argv[3]
    include_d = int(sys.argv[4])
    instance_num = int(sys.argv[5])
    seed(instance_num)

    d_percent = int(percentage)
    
    with open(cube_filename, "r") as cfile:
        data = json.load(cfile)
        cube = data['{}'.format(percentage)]

    cube = cube.split(" ")
    cube_clauses = []
    for val in cube:
        if val != "0":
            cube_clauses.append("\n{} 0".format(val))

    encoding = sys.stdin.readlines()
    clauses = []
    if include_d:
        from math import floor, sqrt
        p_list, q_list, d_list, vars, N = fetch_data(encoding)
        n = len(d_list)

        two_p_list = p_list.copy()
        two_p_list.append(-1)
        two_p_list.reverse()

        two_q_list = q_list.copy()
        two_q_list.append(-1)
        two_q_list.reverse()
        
        two_d_list = d_list.copy()
        two_d_list.append(-1)
        two_d_list.reverse()

        d_list_copy = d_list.copy()
        
        d_list.reverse()

        vars_to_be_added = [two_p_list, two_q_list, two_d_list, d_list]
        addition_list = []
        while len(two_d_list) > 0:
            temp_list = []
            for l in vars_to_be_added:
                if len(l) > 0:
                    v = l.pop(0)
                    if v != -1:
                        temp_list.append(v)
            addition_list.append(temp_list)

        count = 0
        flag = False
        while not flag:
            clauses, vars, flag = add_vars(clauses, vars, addition_list)

        addition_list.reverse()

        sum_val = 2*N + 3
        sum_val = convert_to_binary(sum_val, len(addition_list))

        for i in range(len(addition_list)):
            if sum_val[i] == "1":
                clauses.append("\n{} 0".format(addition_list[i][0]))
            else:
                clauses.append("\n-{} 0".format(addition_list[i][0]))

        with open(data_filename) as fp:
            f = json.load(fp)
            p = f['p']
            q = f['q']
        
        phi = (p-1)*(q-1)
        e = 3
        d_real_val = pow(e, -1, phi) # Needs Python 3.8+
        d_real_bin = convert_to_binary(d_real_val, n)

        set_known_d_bits_count = int(len(d_list_copy) * d_percent/100)
        random_d_bits = sample(d_list_copy, set_known_d_bits_count)

        for bit in random_d_bits:
            if d_real_bin[d_list_copy.index(bit)] == "1":
                clauses.append("\n{} 0".format(bit))
            else:
                clauses.append("\n-{} 0".format(bit))

        dbar_val1 = (2*N + 3) // 3
        dbar_val2 = dbar_val1 - floor(sqrt(2*N))
        dbar_val1 = convert_to_binary(dbar_val1, n)
        dbar_val2 = convert_to_binary(dbar_val2, n)

        dbar_val = ""
        for i in range(n):
            if dbar_val1[i] == dbar_val2[i]:
                dbar_val += dbar_val1[i]
            else:
                break

        for i in range(len(dbar_val)):
            if dbar_val[i] == "1":
                clauses.append("\n{} 0".format(d_list_copy[i]))
            else:
                clauses.append("\n-{} 0".format(d_list_copy[i]))

        l = encoding[14]
        l = l.split(" ")
        l[2] = str(vars)
        l[3] = str(int(l[3]) + len(clauses))
        encoding[14] = " ".join(l) + "\n"

    for line in encoding:
        if line[0]=="p":
            line = line.split(" ")
            line[3] = str(int(line[3]) + len(cube_clauses))
            line = " ".join(line) + "\n"
        print(line, end='')
    for cube in cube_clauses:
        print(cube, end='')
    if clauses:
        for clause in clauses:
            print(clause, end='')
