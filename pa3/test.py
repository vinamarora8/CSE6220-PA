import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('nprocs', type=int, help="Number of CPUs")
parser.add_argument('--n', type=int, default=None, help="Problem size (n)")
parser.add_argument('--debug', type=bool, default=False, help="Debug output")
parser.add_argument('--mul', type=bool, default=False, help="Will test mat-vec-mul if true")
args = parser.parse_args()

nprocs = args.nprocs
n = args.n
assert(int(np.sqrt(nprocs))**2 == nprocs)

if n is None:
    n = np.random.randint(4, 100)
A = np.random.random((n, n)) * 10
# Make the matrix diagonally dominant
for i in range(n):
    diagonal_element = A[i, i]
    row_sum = np.sum(np.abs(A[i, :])) - np.abs(diagonal_element)
    if diagonal_element <= row_sum:
        A[i, i] += row_sum - diagonal_element + 1

x = np.random.random(n) * 100

inp_mat_fname = 'inp_mat.txt'
inp_vec_fname = 'inp_vec.txt'
out_fname = 'out_vec.txt'

if args.mul:
    y = np.dot(A, x)
else:
    y = np.linalg.solve(A, x)

def write_matrix(A, fname):
    n = A.shape[0]
    with open(fname, 'w') as f:
        f.write(f'{n}\n')
        for i in range(n):
            for j in range(n):
                f.write(f'{A[i, j]} ')
            f.write('\n')

def write_vector(x, fname):
    with open(fname, 'w') as f:
        for i in range(len(x)):
            f.write(f'{x[i]} ')
        f.write('\n')

def read_vector(fname, n):
    with open(fname, 'r') as f:
        text = f.readlines()[0].strip()
    text_double = [float(x) for x in text.split(' ')]
    ans = np.array(text_double)
    return ans

def run_program():
    cmd = f'mpirun -np {nprocs} --oversubscribe src/pjacobi {inp_mat_fname} {inp_vec_fname} {out_fname}'
    print(cmd)
    os.system(cmd)

write_matrix(A, inp_mat_fname)
write_vector(x, inp_vec_fname)
run_program()


ans = read_vector('out_vec.txt', n)
err = y - ans
err_norm = np.linalg.norm(err)
criteria = 1e-9 * n

if args.debug:
    print('Expected')
    print(y)

    print('Actual')
    print(ans)

    print(f'Error:')
    print(err)

print(f'Error norm: {err_norm}')
print(f'Criteria: {criteria}')
if (err_norm < criteria):
    print('PASS')
    exit(0)
else:
    print('FAIL')
    exit(1)
