import numpy as np
import os
import matplotlib.pyplot as plt
import subprocess


start = 7
end = 20
steps = end-start+1
n = np.logspace(start, end, num=steps, base=2).astype(int)

processors = [2, 4, 8, 12, 24]

for p in processors:
    exe_times = []
    for m in n:

        cmd_base = 'python test.py {} {}'.format(p, int(m))
        print(cmd_base)
        status, output = subprocess.getstatusoutput(cmd_base)
        
        lines = output.splitlines()
        print(float(lines[0]))
        exe_times.append(float(lines[0]))
    
    print(n)        
    print(exe_times)

    plt.clf()
    plt.xscale('log', base=2)
    plt.plot(n, exe_times)
    plt.xlabel("problem size")
    plt.ylabel("execution time in ms")
    plt.title("Performance analysis for p = {}".format(p))

    # Save the plot to a file named "my_plot.png"
    plt.savefig("plots/perf_{}.png".format(p))


