import numpy as np
import os
import subprocess


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

def run_program(p, inp_mat_fname, inp_vec_fname, out_fname):

    cmd = f'mpirun -np {p} src/pjacobi {inp_mat_fname} {inp_vec_fname} {out_fname}'
    try:
        output = subprocess.check_output(cmd, shell=True)
        output = str(output).strip().split()
        runtime = float(output[1])
    except subprocess.CalledProcessError as e:
        print(e.output)

    print(f'Runtime: {runtime}')
    return runtime


def one(n, p, debug=False, mul=False, passfail=False):

    assert(int(np.sqrt(p))**2 == p)

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

    if mul:
        y = np.matmul(A, x)
    else:
        y = np.linalg.solve(A, x)

    inp_mat_fname = 'inp_mat.txt'
    inp_vec_fname = 'inp_vec.txt'
    out_fname = 'out_vec.txt'

    write_matrix(A, inp_mat_fname)
    write_vector(x, inp_vec_fname)
    runtime = run_program(p, inp_mat_fname, inp_vec_fname, out_fname)

    ans = read_vector('out_vec.txt', n)
    err = y - ans
    err_norm = np.linalg.norm(err)
    criteria = 1e-9 * n

    if debug:
        print('Expected')
        print(y)

        print('Actual')
        print(ans)

        print(f'Error:')
        print(err)

    if passfail:
        print(f'Error norm: {err_norm}')
        print(f'Criteria: {criteria}')
        if (err_norm < criteria):
            print('PASS')
            exit(0)
        else:
            print('FAIL')
            exit(1)

    assert(err_norm < criteria)

    return runtime


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('nprocs', type=int, help="Number of CPUs")
    parser.add_argument('--n', type=int, default=None, help="Problem size (n)")
    parser.add_argument('--debug', type=bool, default=False, help="Debug output")
    parser.add_argument('--mul', type=bool, default=False, help="Will test mat-vec-mul if true")
    args = parser.parse_args()
    nprocs = args.nprocs

    one(args.n, args.nprocs, debug=args.debug, mul=args.mul, passfail=True)
