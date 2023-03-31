import numpy as np
import os
import matplotlib.pyplot as plt
import subprocess

p_max = 4
start = 7
end = 20
steps = end-start+1
n = np.logspace(start, end, num=steps, base=2).astype(int)
x = range(1, p_max+1)

processors = range(1, p_max+1)

for m in n:
    exe_times = []
    for p in processors:
    
        cmd_base = 'python test.py {} {}'.format(p, int(m))
        print(cmd_base)
        status, output = subprocess.getstatusoutput(cmd_base)
        
        lines = output.splitlines()
        print(float(lines[0]))
        exe_times.append(float(lines[0]))
    
    
    print(exe_times)

    plt.clf()
    plt.plot(x, exe_times)
    plt.xlabel("num of processors")
    plt.ylabel("execution time in ms")
    plt.title("Performance analysis for m = {}".format(int(m)))

    # Save the plot to a file named "my_plot.png"
    plt.savefig("plots/perf_{}.png".format(int(np.log2(m))))


