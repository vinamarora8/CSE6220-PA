#!/bin/bash

set -e

make

for np in {1..24}; do
    echo -n "$np, "
    mpiexec -n $np ./int_calc 1000000
done
