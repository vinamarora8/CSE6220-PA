
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import io 
import sys
import os
import re

from test import one



'''
get runtime vs p
for n = 10, 100, 1000, 10000, 100000, 1000000
for p = 1, 4, 9, 16, 25, 36, 49, 64
plot runtime vs p for each n
'''
def runtime_vs_p():
    procs = [1, 4, 9, 16, 25, 36, 49, 64]
    nums = [10, 100, 1000, 10000, 100000, 1000000]
    avg_runtimes = []
    for n in nums:
        temp_runtimes = []
        print('Starting for n = ', n)
        for p in procs:
            sum = 0
            # find the mean of 5 runs
            for i in range(5):
                # run test.py with n, p
                cmd = f'python test.py {p} --n {n}'
                output = subprocess.check_output(cmd, shell=True)
                output_lines = output.decode('utf-8').split('\n')
                pattern = r"Runtime: (\d+\.\d+)"
                match = re.search(pattern, output_lines[0])
                sum += float(match.group(1))
            temp_runtimes.append(sum/5)
        print('Done for n = ', n)
        avg_runtimes.append(temp_runtimes)

    # plot runtime vs p for each n 
    # juxtapose all plots
    for i in range(len(nums)):
        plt.plot(procs, avg_runtimes[i], label=f"n = {nums[i]}")
        plt.scatter(procs, avg_runtimes[i])
    plt.grid()
    plt.yscale('log')
    plt.xticks(procs, procs)
    plt.xlabel('Number of processes (p)')
    plt.ylabel('Run-time (ms)')
    plt.legend()
    plt.savefig('runtime_vs_p.png')
    plt.clf()


'''
get runtime vs n
for n = 10, 100, 1000, 10000, 100000, 1000000
for p = 1, 4, 9, 16, 25, 36, 49, 64
plot runtime vs n for each p
'''
def runtime_vs_n():

    procs = [1, 4, 9, 16, 25, 36, 49, 64]
    nums = [10 , 100, 1000, 10000, 100000, 1000000]
    avg_runtimes = []
    for p in procs:
        temp_runtimes = []
        print('Starting for p = ', p)
        for n in nums:
            sum = 0
            # find the mean of 5 runs
            for i in range(5):
                # run test.py with n, p
                cmd = f'python test.py {p} --n {n}'
                output = subprocess.check_output(cmd, shell=True)
                output_lines = output.decode('utf-8').split('\n')
                pattern = r"Runtime: (\d+\.\d+)"
                match = re.search(pattern, output_lines[0])
                sum += float(match.group(1))
            temp_runtimes.append(sum/5)
        print('Done for p = ', p)
        avg_runtimes.append(temp_runtimes)

    # plot runtime vs n for each p
    # juxtapose all plots
    for i in range(len(procs)):
        plt.plot(nums, avg_runtimes[i], label=f"p = {procs[i]}")
        plt.scatter(nums, avg_runtimes[i])
    plt.grid()
    plt.yscale('log')
    plt.xscale('log', base=10)
    plt.xlabel('Number of elements (n)')
    plt.ylabel('Run-time (ms)')
    plt.legend()
    plt.savefig('runtime_vs_n.png')
    plt.clf()

print('Running runtime_vs_n()')
runtime_vs_n()
print('Running runtime_vs_p()')
runtime_vs_p()