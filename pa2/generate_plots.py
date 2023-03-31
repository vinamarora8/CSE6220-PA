import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import random 
import time

script_path = sys.argv[1]

def get_times():

    # problem size 
    

    n = np.logspace(32, 1048576, 2)

    for m in n:
        # generate wrt to power of 2

        processors = np.logspace(1, (m/2), 2)
        
        exe_times = []
        for p in processors:    
            
            random_integers = [random.randint(1, 10*m) for _ in range(m)]

            # Write the integers to a file
            with open("input.txt", "w") as file:
                file.write(str(m) + "\n")  
                file.write(" ".join(str(x) for x in random_integers)) 
        
            cmd_base = 'mpiexec -np {} ./pqsort  input.txt output.txt'.format(p)

            print(cmd_base)

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

