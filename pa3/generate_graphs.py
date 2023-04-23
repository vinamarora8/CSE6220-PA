
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import re
from test import one


procs = [4, 9, 16, 25, 36, 49]
nums = [10, 50, 100, 300, 500, 1000] 


def util(n, p):
    # run test.py with n, p  
    cmd = f'python test.py {p} --n {n}'
    print(cmd)
    
    # test once so that the output is not messed up
    try:
        output = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(e.output)
    
    count = 0
    sum = 0
    temp_runtime = 0

    for i in range(3):
        output = subprocess.check_output(cmd, shell=True)
        output_lines = output.decode('utf-8').split('\n')
        pattern = r"Runtime: (\d+\.\d+)"
        match = re.search(pattern, output_lines[0])
        if match != None:
            count += 1
            sum += float(match.group(1))
    if count != 0:
        temp_runtime = sum/count
    else:
        print("For n = ", n, " and p = ", p, " count = ", count)
    return temp_runtime


def runtime_vs_p():

    avg_runtimes = []
    for n in nums:
        temp_runtimes = []
        for p in procs:
            temp_runtimes.append(util(n, p))
        print(temp_runtimes)
    
    avg_runtimes.append(temp_runtimes)
    print(avg_runtimes)

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


def runtime_vs_n():
    avg_runtimes = []
    for p in procs:
        temp_runtimes = []
        for n in nums:
            temp_runtimes.append(util(n, p))
        print(temp_runtimes)
    
    avg_runtimes.append(temp_runtimes)
    print(avg_runtimes)

    # plot runtime vs n for each p
    # juxtapose all plots
    for i in range(len(procs)):
        plt.plot(nums, avg_runtimes[i], label=f"p = {procs[i]}")
        plt.scatter(nums, avg_runtimes[i])
    plt.grid()
    plt.yscale('log')
    plt.xticks(nums, nums)
    plt.xlabel('Number of elements (n)')
    plt.ylabel('Run-time (ms)')
    plt.legend()
    plt.savefig('runtime_vs_n.png')
    plt.clf()

print('Running runtime_vs_n()')
runtime_vs_n()
print('Running runtime_vs_p()')
runtime_vs_p()