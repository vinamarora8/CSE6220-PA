import sys
import os
import random
import numpy as np
import argparse

def gen_input(size, op_fname, log=False):
    int_min = -2**31
    int_max = (2**31)-1
    nums = [random.randint(int_min, int_max) for _ in range(size)] 

    with open(op_fname, 'w') as f:
        f.write(f'{size}\n')
        for n in nums:
            f.write(f'{n} ')
        f.write('\n')

    if log:
        print(f'inp : {nums}')
    return nums


def check_output(fname, log=False):
    with open(fname) as f:
        text = f.readlines()

    time = float(text[1].strip())
    text = text[0].strip().split(' ')
    nums = [int(x) for x in text]

    arr = np.array(nums)
    diff = arr[1:] - arr[:-1]
    check = (diff >= 0).all()
    assert check

    if log:
        print(f'out : {nums}')
        print(f'diff: {diff}')

    return time


def run_prog(ip_fname, op_fname, num_procs):
    cmd = f'mpirun -np {num_procs} pqsort {ip_fname} {op_fname}'
    os.system(cmd)


parser = argparse.ArgumentParser()
parser.add_argument('num_procs', type=int)
parser.add_argument('size', type=int)
parser.add_argument('--log', action='store_true')
args = parser.parse_args()

size = args.size
num_procs = args.num_procs
log = args.log
ip_fname = "i.txt"
op_fname = "o.txt"


gen_input(size, ip_fname, log)
run_prog(ip_fname, op_fname, num_procs)
time = check_output(op_fname, log)
print(f'{time}')
