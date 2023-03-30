How it works:
-------------
The summation of n-terms is split p summations of (n/p) terms.
Rank 0 reads the value of `n` from command line args and broadcasts
it to all ranks.

Each process then computes its own iteration start and end points,
and performs the summation in a `local_sum` variable. Finally, all
individual local sums are reduced by the MPI_SUM operation on rank 0.


Machine:
--------
coc-ice-multi cluster, with resource string "nodes=128:ppn=1"
