import subprocess
import matplotlib.pyplot as plt
import numpy as np
from test import gen_input, check_output, run_prog

ip_fname = "i.txt"
op_fname = "o.txt"

def get_avg_time(n, p, runs):
    times = []
    for i in range(runs):
        gen_input(n, ip_fname)
        run_prog(ip_fname, op_fname, p)
        times.append(check_output(op_fname))
    avg_time = np.mean(times)
    std      = np.std(times)
    return avg_time, std

n = [100, 1000, 10000, 100000, 1000000, 10000000]
p = 8

exe_times = []
exe_times_std = []
for n_i in n:
    
    time_avg, time_std = get_avg_time(n_i, p, 8)
    exe_times.append(time_avg)
    exe_times_std.append(time_std)
    print(f'{n_i}, {p}, {time_avg}, {time_std}')

print(exe_times)
print(exe_times_std)

for i in range(len(n)):
    print(f'{n[i]}, {p}, {exe_times[i]}, {exe_times_std[i]}')

'''
plt.clf()
plt.xscale('log', base=2)
plt.xlabel("problem size")
plt.ylabel("execution time in ms")
plt.title("Performance analysis for p = {}".format(p))

# Save the plot to a file named "my_plot.png"
plt.savefig("plots/perf_{}.png".format(n))
'''