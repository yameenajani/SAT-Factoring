import sys
from ast import literal_eval
import linecache
import math
from sage.all import *
import json
from random import randint, sample

# fetch_lists reads the instance and fetches the lists of variables corresponding to p, q and d
def fetch_data(input_filename):
    q_list = linecache.getline(input_filename, 12)
    q_list = q_list.split(": ")
    q_list = literal_eval(q_list[1])

    p_list = linecache.getline(input_filename, 11)
    p_list = p_list.split(": ")
    p_list = literal_eval(p_list[1])
    make_msb_0_len = len(p_list) - len(q_list)
    p_list = p_list[make_msb_0_len:]

    N_size = linecache.getline(input_filename, 9)
    N_size = N_size.split(" ")
    N = int(N_size[5])
    N_size = literal_eval(N_size[6])
    N_size = len(N_size)

    num_variables = linecache.getline(input_filename, 15)
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

def convert_to_binary(dec_val):
    return bin(dec_val).split("b")[1]

def main():
    bitsize = sys.argv[1]
    percentage = sys.argv[2]
    instance_number = sys.argv[3]
    
    d_percent = int(percentage)

    data_filename = "../data/{}/data_{}.json".format(bitsize, instance_number)
    instance_filename = "../instances/{}/temp_withd/{}/instance_{}.cnf".format(bitsize, percentage, instance_number)

    p_list, q_list, d_list, vars, N = fetch_data(instance_filename)
    
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

    clauses = []

    count = 0
    flag = False
    while not flag:
        clauses, vars, flag = add_vars(clauses, vars, addition_list)

    addition_list.reverse()

    sum_val = 2*N + 3
    sum_val = convert_to_binary(sum_val)
    pad_len = len(addition_list) - len(sum_val)
    for _ in range(pad_len):
        sum_val = "0" + sum_val

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
    e = ZZ(3)
    d_real_val = e.inverse_mod(phi)
    d_real_bin = convert_to_binary(d_real_val)
    pad_len = len(d_list_copy) - len(d_real_bin)
    for _ in range(pad_len):
        d_real_bin = "0" + d_real_bin

    set_known_d_bits_count = int(len(d_list_copy) * d_percent/100)
    random_d_bits = sample(d_list_copy, set_known_d_bits_count)

    for bit in random_d_bits:
        if d_real_bin[d_list_copy.index(bit)] == "1":
            clauses.append("\n{} 0".format(bit))
        else:
            clauses.append("\n-{} 0".format(bit))

    dbar_val1 = (2*N + 3) // 3
    dbar_val2 = dbar_val1 - 2*floor(sqrt(N))
    dbar_val1 = convert_to_binary(dbar_val1)
    dbar_val2 = convert_to_binary(dbar_val2)
    dbar_val = ""
    for i in range(min(len(dbar_val1), len(dbar_val2))):
        if dbar_val1[i] == dbar_val2[i]:
            dbar_val += dbar_val1[i]
        else:
            break

    for _ in range(pad_len):
        dbar_val = "0" + dbar_val

    for i in range(len(dbar_val)):
        if dbar_val[i] == "1":
            clauses.append("\n{} 0".format(d_list_copy[i]))
        else:
            clauses.append("\n-{} 0".format(d_list_copy[i]))

    with open(instance_filename, "r") as f:
        lines = f.readlines()

    l = lines[14]
    l = l.split(" ")
    l[2] = str(vars)
    l[3] = str(int(l[3]) + len(clauses))
    lines[14] = " ".join(l) + "\n"

    with open(instance_filename, "w") as f:
        for line in lines:
            f.write(line)
        for clause in clauses:
            f.write(clause)


if __name__ == '__main__':
    main()