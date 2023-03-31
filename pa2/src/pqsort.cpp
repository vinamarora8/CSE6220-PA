#include <iostream>
#include <fstream>
#include <mpi.h>

#define ROOT 0
//#define DEBUG_BLOCK_DIST

void serial_sort(int *inp, int low, int high);
void parallel_qsort(int *inp, int len, int global_len, int seed, MPI_Comm comm);
int distribute_input(const char *fname, int *&local_inp, MPI_Comm comm);
void gather_output(int *local_arr, int local_len, std::ofstream &fstream, MPI_Comm comm);

int main(int argc, char *argv[])
{

    char *in_fname = argv[1];
    char *op_fname = argv[2];

    int p, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Rank 0 reads input file and block distributes
    int local_len, *local_inp;
    local_len = distribute_input(in_fname, local_inp, MPI_COMM_WORLD);

    // Timing start
    double starttime = MPI_Wtime();

    //parallel_qsort(local_inp, local_len, 0, MPI_COMM_WORLD);

    // Timing end
    double runtime = (MPI_Wtime() - starttime) * 1000.0;

    // Print to output file
    std::ofstream opfile;
    if (rank == ROOT) opfile.open(argv[2]);
    gather_output(local_inp, local_len, opfile, MPI_COMM_WORLD);
    if (rank == ROOT) opfile << std::fixed << std::setprecision(6) << runtime << std::endl;
    if (rank == ROOT) opfile.close();

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

void parallel_qsort(int *inp, int len, int global_len, int seed, MPI_Comm comm)
{

    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    if (p == 1)
    {
        serial_sort(inp, 0, len - 1);
        return;
    }

    // Choose pivot (random, with same seed) and broadcast
    // TODO 
    int pivot_index = -1;
    int pivot = -1;

    std::srand(seed);
    pivot_index = std::rand() % global_len;

    // find which processor has the pivot
    int pivot_rank = -1;
    if ( (pivot_index >= rank * global_len / p ) && (pivot_index < (rank + 1) * global_len / p) ) {
        pivot = inp[pivot_index];
    }
    MPI_Bcast(&pivot, 1, MPI_INT, rank, comm);


    // Local partition, and count sizes
    // TODO
    int local_low_len = 0, local_high_len = 0;

    int i = - 1;
    for (int j = 0; j <= len-1; j++) { 
        if (inp[j] <= pivot) {
            i++;
            std::swap(inp[i], inp[j]);
        }
    }

    local_low_len = i;
    local_high_len = len - i;

    // All-gather on lengths
    // TODO
    int* low_len = new int[p];
    int* high_len = new int[p];
    int new_global_len;

    MPI_Allgather(&local_low_len, 1, MPI_INT, low_len, 1, MPI_INT, comm);
    MPI_Allgather(&local_high_len, 1, MPI_INT, high_len, 1, MPI_INT, comm);

    // Compute low/high partitions and how many procs to assign for each
    // Create new communicators
    // TODO
    int p_low, p_high;
    MPI_Comm comm_low, comm_high, my_comm;

    // Compute communication indices + All-to-all
    // TODO
    int *new_arr, new_len;

    // Compute new seeds
    int new_seed = std::rand();

    if (new_len > 0)
        parallel_qsort(new_arr, new_len, new_global_len, new_seed, my_comm);
}


/**
 * Read array from input file on rank 0 and block-distribute to all processes 
 * in communicator.
 * Each process gets either ceil(n/p) or floor(n/p) values.
 * 
 * @return local_inp (in-place)
 * @return length of local_inp
 */
int distribute_input(const char *fname, int *&local_inp, MPI_Comm comm)
{
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    int len, *inp;
    if (rank == 0)
    {
        // Load input file into array
        std::ifstream file(fname);

        file >> len; // first line is length of array
        // remaining lines are the array separated by spaces
        inp = (int *) malloc(len * sizeof(int));
        for (int i = 0; i < len; i++)
            file >> inp[i];
        file.close();

#ifdef DEBUG_INP_READ
        std::cout << len << std::endl;
        for (int i = 0; i < len; i++)
            std::cout << inp[i] << " ";
        std::cout << std::endl;
#endif
    }


    // Decide split counts
    int counts[p], displs[p];
    if (rank == 0)
    {
        for (int i = 0; i < p; i++)
        {
            counts[i] = (len / p) + (i < (len % p));
            displs[i] = i == 0 ? 0 : displs[i-1] + counts[i-1];
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


/**
 * Gather final arrays from all processors and write to fstream
 */
void gather_output(int *local_arr, int local_len, std::ofstream &fstream, MPI_Comm comm)
{
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    int total_len;
    int counts[p], displs[p];

    // Gather total length, counts, and compute displacements
    MPI_Reduce(&local_len, &total_len, 1, MPI_INT, MPI_SUM, ROOT, comm);
    MPI_Gather(
        &local_len, 1, MPI_INT,
        (void *) counts, 1, MPI_INT,
        ROOT, comm
    );
    if (rank == ROOT)
        for (int i = 0; i < p; i++)
            displs[i] = i == 0 ? 0 : displs[i-1] + counts[i-1];

    // Gather outputs
    int *comb_array = (int *) malloc(total_len * sizeof(int));
    MPI_Gatherv(
        local_arr, local_len, MPI_INT,
        comb_array, counts, displs, MPI_INT,
        ROOT, comm
    );

    if (rank == ROOT)
    {
        for (int i = 0; i < total_len; i++)
        {
            std::cout << comb_array[i] << " ";
            fstream << comb_array[i] << " ";
        }
        std::cout << std::endl;
        fstream << std::endl;
    }

    free(comb_array);
}

