import subprocess
import matplotlib.pyplot as plt


n=10000
p = range(1, 25)


exe_times = []
for p_i in p:
    
    cmd_base = 'python test.py {} {}'.format(p_i, n)
    print(cmd_base)
    status, output = subprocess.getstatusoutput(cmd_base)
    
    lines = output.splitlines()
    print(float(lines[0]))
    exe_times.append(float(lines[0]))
          
print(exe_times)

plt.clf()
plt.xscale('log', base=2)
plt.plot(p, exe_times)
plt.xlabel("problem size")
plt.ylabel("execution time in ms")
plt.title("Performance analysis for p = {}".format(p))

# Save the plot to a file named "my_plot.png"
plt.savefig("plots/perf_{}.png".format(n))


