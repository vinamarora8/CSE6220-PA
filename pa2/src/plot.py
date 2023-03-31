import numpy as np
import os
import matplotlib.pyplot as plt
import subprocess


n = np.logspace(10, 20, num=8, base=2)
x = range(1, 9)
for m in n:
    exe_times = []
    for p in range(1, 9):

        cmd_base = 'python test.py {} {}'.format(p, int(m))
        print(cmd_base)
        status, output = subprocess.getstatusoutput(cmd_base)
        

        lines = output.splitlines()
        print(lines[0])
        exe_times.append(lines[0])

    plt.clf()
    plt.plot(x, exe_times)
    plt.xlabel("number of processors")
    plt.ylabel("execution time")
    plt.title("Performance analysis")

    # Save the plot to a file named "my_plot.png"
    plt.savefig("perf_{}.png".format(m))


