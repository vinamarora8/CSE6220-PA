#include <iostream>
#include <fstream>
#include <mpi.h>

void serial_sort(int *inp, int low, int high);
void parallel_qsort(int *inp, int len, int seed, MPI_Comm comm);

int main(int argc, char *argv[])
{
    char *inp_fname = argv[1];
    char *out_fname = argv[2];

    int p, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Rank 0 reads input file and block distributes
    int local_len, *local_inp;

    // Timing start
    double starttime = MPI_Wtime();

    parallel_qsort(local_inp, local_len, 0, MPI_COMM_WORLD);

    // Timing end
    double runtime = MPI_Wtime() - starttime;

    // Print to output file
    // TODO

    MPI_Finalize();
    return 0;
}

// A helper function to swap two elements
void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// A helper function to partition the array around a pivot
int partition(int *arr, int low, int high) {
    int pivot = arr[high];
    int i = low - 1;
    for (int j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

// The serial sort function (Quick sort implemented)
void serial_sort(int *arr, int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);
        serial_sort(arr, low, pi - 1);
        serial_sort(arr, pi + 1, high);
    }
}

void parallel_qsort(int *inp, int len, int seed, MPI_Comm comm)
{

    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    if (p == 1)
    {
        serial_sort(inp, 0, len - 1);
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
    int *new_low, *new_high;

    // Compute new seeds
    int seed_low, seed_high;

    parallel_qsort(new_low, new_low_len, seed_low, comm_low);
    parallel_qsort(new_high, new_high_len, seed_high, comm_high);
}
