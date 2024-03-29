import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import random 
import time

def get_times():

    # problem size 
    

    n = np.logspace(3, 10, num=8, base=2).astype(int)

    for m in n:
        
        # generate wrt to power of 2
        num=int(np.log2(int(m)))
        processors = np.logspace(1, num-1, num=(num-1), base=2).astype(int)
        
        exe_times = []
        for p in processors:    
            
            random_integers = [random.randint(1, 10*m) for _ in range(m)]

            # Write the integers to a file
            with open("input.txt", "w") as file:
                file.write(str(m) + "\n")  
                file.write(" ".join(str(x) for x in random_integers)) 
        
            cmd_base = 'mpiexec -np {} ./pqsort  input.txt output.txt'.format(p)

            print(cmd_base)

            # add start, end time calculation in cpp and redirect them to tmp
            start_time = time.time()
            os.system(cmd_base)
            end_time = time.time()
            time_taken = end_time - start_time
            print("Time taken: {:.6f} seconds".format(time_taken))

            exe_times.append(time_taken)

        plt.plot(p, exe_times)
        plt.xlabel("number of processors")
        plt.ylabel("execution time")
        plt.title("Performance analysis")

        # Save the plot to a file named "my_plot.png"
        plt.savefig("perf_{}_{}.png".format(m, p))
    
    return


# mpicxx pqsort.cpp -o pqsort
# mpiexec -np 4 ./pqsort input.txt output.txt

get_times()

