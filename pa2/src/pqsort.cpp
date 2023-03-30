#include <iostream>
#include <fstream>
#include <mpi.h>

#define ROOT 0

void serial_sort(int *inp, int low, int high);
void parallel_qsort(int *inp, int len, int seed, MPI_Comm comm);
int distribute_input(const char *fname, int *local_inp, MPI_Comm comm);


int main(int argc, char *argv[])
{
    int p, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Rank 0 reads input file and block distributes
    int local_len, *local_inp;
    local_len = distribute_input(argv[1], local_inp, MPI_COMM_WORLD);

    // Timing start
    double starttime = MPI_Wtime();

    //parallel_qsort(local_inp, local_len, 0, MPI_COMM_WORLD);

    // Timing end
    double runtime = MPI_Wtime() - starttime;

    // Print to output file
    if (rank == 0)
    {
        char *out_fname = argv[2];
    }

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
    MPI_Comm comm_low, comm_high, my_comm;

    // Compute communication indices + All-to-all
    // TODO
    int *new_arr, new_len;

    // Compute new seeds
    int seed_low, seed_high;

    if (new_len > 0)
        parallel_qsort(new_arr, new_len, seed, my_comm);
}


/**
 * Read array from input file on rank 0 and block-distribute to all processes 
 * in communicator.
 * Each process gets either ceil(n/p) or floor(n/p) values.
 * 
 * @return local_inp (in-place)
 * @return length of local_inp
 */
int distribute_input(const char *fname, int *local_inp, MPI_Comm comm)
{
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    int len, *inp;
    if (rank == 0)
    {
        // Load input file into array
        std::ifstream file(fname);

        file >> len;
        inp = (int *) malloc(len * sizeof(int));
        for (int i = 0; i < len; i++)
            file >> inp[i];

#ifdef DEBUG_INP_READ
        std::cout << len << std::endl;
        for (int i = 0; i < len; i++)
            std::cout << inp[i] << " ";
        std::cout << std::endl;
#endif
    }


    // Decide split counts
    int counts[p], displs[p];
    displs[0] = 0;
    if (rank == 0)
    {
        for (int i = 0; i < p; i++)
        {
            counts[i] = (len / p) + (i < (len % p));
            displs[i+1] = displs[i] + counts[i];
        }
    }


    // Communicate lengths, create buffers, and copy data
    int local_len;
    MPI_Scatter(
        (int *) counts, 1, MPI_INT, 
        &local_len, 1, MPI_INT, 
        ROOT, comm
    );

    local_inp = (int *) malloc(local_len * sizeof(int));

    MPI_Scatterv(
        (void *) inp, (int *) counts, (int *) displs, MPI_INT,
        (void *) local_inp, local_len, MPI_INT,
        ROOT, comm
    );

    if (rank == 0)
        free(inp);

#ifdef DEBUG_BLOCK_DIST
    for (int i = 0; i < p; i++)
    {
        if (i == rank)
        {
            std::cout << "Rank " << i << ": ";
            for (int j = 0; j < local_len; j++)
                std::cout << local_inp[j] << " ";
            std::cout << std::endl;
        }

        MPI_Barrier(comm);
    }
#endif
    
    return local_len;
}
