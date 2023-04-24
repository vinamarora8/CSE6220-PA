import numpy as np
import matplotlib.pyplot as plt

data_fname = "perf.txt"
data = np.loadtxt(data_fname)

# remove data points with p not being a power of 2
data = [x for x in data if (2**(np.log2(x[1])) == x[1] and x[0] >= 16)]
data = np.array(data)

n = np.unique(data[:, 0]).astype(int)
p = np.unique(data[:, 1]).astype(int)


# runtimes vs n
for i in range(len(p)):
    
    # if p[i] is power of 2
    times_n = [x for x in data if x[1] == p[i]]
    
    ns = [x[0] for x in times_n]
    times = [x[2]/1000. for x in times_n]

    plt.plot(ns, times, label=f"p = {p[i]}")
    plt.scatter(ns, times, s=8)

plt.legend()
plt.xlabel('Linear size of matrix A (n)')
plt.ylabel('Run-time (s)')
plt.grid()
plt.yscale('log')
plt.xscale('log', base=2)
plt.xticks(n)
plt.show()

# runtimes vs p
for i in range(len(n)):
    if n[i] < 64 or n[i] > 1000:
        continue
    times_p = [x for x in data if x[0] == n[i]]
    ps = [x[1] for x in times_p]
    times = [x[2]/1000. for x in times_p]

    plt.plot(ps, times, label=f"n = {n[i]}")
    plt.scatter(ps, times, s=8)

    p_ideal = np.argmin(times)
    plt.scatter(ps[p_ideal], times[p_ideal], s=16, c='k')

plt.legend()
plt.xlabel('Number of processes (p)')
plt.ylabel('Run-time (s)')
plt.grid()
plt.yscale('log')
plt.xscale('log', base=2)
plt.xticks(p)
plt.show()
