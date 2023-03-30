#include <iostream>
#include <fstream>
#include <mpi.h>

int serial_sort(int *inp, int len);
int parallel_qsort(int *inp, int len, MPI_Comm comm);

int main(int argc, char *argv[])
{
    char *inp_fname = argv[1];
    char *out_fname = argv[2];

    int p, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Rank 0 reads input file and block distributes
    int local_len, *local_inp;

    // Timing start
    // TODO

    parallel_qsort(local_inp, local_len, 0, MPI_COMM_WORLD);

    // Timing end
    // TODO

    // Print to output file
    // TODO
}


int serial_sort(int *inp, int len)
{
    // TODO
}


int parallel_qsort(int *inp, int len, int seed, MPI_Comm comm)
{

    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    if (p == 1)
    {
        serial_sort(inp, len);
    }

    // Choose pivot (random, with same seed) and broadcast
    // TODO
    int pivot;

    // Local partition, and count sizes
    // TODO
    int *local_low, *local_high;
    int local_low_len, local_high_len;

    // All-gather on lengths
    // TODO
    int low_len[p], high_len[p];

    // Compute low/high partitions and how many procs to assign for each
    // Create new communicators
    // TODO
    int p_low, p_high;
    MPI_Comm comm_low, comm_high;

    // Compute communication indices + All-to-all
    // TODO
    int new_low_len, new_high_len;
    int *new_lowx, *new_high;

    // Compute new seeds
    int seed_low, seed_high;

    parallel_qsort(new_low, new_low_len, seed_low, comm_low);
    parallel_qsort(new_high, new_high_len, seed_high, comm_high);
}
