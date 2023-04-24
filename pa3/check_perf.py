from test import one
import numpy as np

p = [4, 9, 16, 25, 36, 49, 64]
n = 2**np.arange(2, 10)

runtimes = np.zeros((len(p), len(n)))

op_fname = "perf.txt"

def avg_runtime(n, p, count):
    times = []
    for i in range(count):
        times.append(one(n, p))
    avg_time = np.mean(times)
    print(f'Average @n={n}, p={p}: {avg_time}')
    return avg_time

for i in range(len(p)):
    for j in range(len(n)):
        runtimes[i, j] = avg_runtime(n[j], p[i], 4)
        with open(op_fname, "a") as f:
            f.write(f'{n[j]} {p[i]} {runtimes[i, j]}\n')
