import numpy as np
import matplotlib.pyplot as plt

avg_times = [
    [147.6327629, 69.19775275, 32.79459413, 18.15306975, 12.055775, 7.145247, 5.9651225],
    [12.36259175, 10.01818363, 5.817194, 2.128753375, 1.443805875, 1.934684125, 2.812175625],
    [1.444793, 1.362932625, 0.509955625, 0.73752075, 0.99344025, 1.64490475, 2.758881],
    [0.119682375, 0.145497, 0.309722625, 0.50910925, 0.761663, 1.4186195, 2.765583375],
]

std_times = [
    [3.758634528, 43.89334799, 26.81701362, 8.934889785, 3.432732428, 2.744981212, 2.058209039],
    [0.172541211, 1.290836792, 3.656523285, 0.809711772, 0.46915301, 0.592228346, 0.776353103],
    [0.33186929, 0.538924744, 0.331303648, 0.152268943, 0.208262388, 0.28151883, 0.303721581],
    [0.03918212, 0.045688978, 0.060613612, 0.125246768, 0.246744454, 0.179781335, 0.445406179],
]

m = [1000000, 100000, 10000, 1000]
p = [1, 2, 4, 8, 16, 32, 64]

for i in range(len(m)):
    #plt.errorbar(p, avg_times[i], std_times[i], capsize=5.0)
    plt.plot(p, avg_times[i])
    plt.scatter(p, avg_times[i], label=f"n = {m[i]}")

plt.grid()
plt.yscale('log')
plt.xscale('log', base=2)
plt.xlabel('Number of processors')
plt.ylabel('Run-time (ms)')
plt.legend()
plt.show()


for i in range(0,len(p),2):
    plt.plot(m, [avg_times[x][i] for x in range(len(m))], label=f"p = {p[i]}")
    plt.scatter(m, [avg_times[x][i] for x in range(len(m))])
plt.grid()
plt.yscale('log')
plt.xscale('log', base=10)
plt.xlabel('Problem size (n)')
plt.ylabel('Run-time (ms)')
plt.legend()
plt.show()