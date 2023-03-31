#include <iostream>
#include <algorithm>
#include <fstream>
#include <mpi.h>

int serial_sort(int *inp, int len);
int parallel_qsort(int *inp, int len, int seed, MPI_Comm comm);

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
    double starttime = MPI_Wtime();

    parallel_qsort(local_inp, local_len, 0, MPI_COMM_WORLD);

    // Timing end
    double runtime = MPI_Wtime() - starttime;

    // Print to output file
    // TODO
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

// Computes the number of overlaps between two ranges.
int count_overlaps(int start1, int end1, int start2, int end2) {
    int max_start = std::max(start1, start2);
    int min_end = std::min(end1, end2);
    int overlap_count = std::max(0, min_end - max_start + 1);
    return overlap_count;
}

int parallel_qsort(int *inp, int len, int seed, MPI_Comm comm)
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
    int local_low_len = 0, local_high_len = 0;

    int i = - 1;
    for (int j = 0; j <= len-1; j++) { 
        if (inp[j] <= pivot) {
            i++;
            local_low_len++;
            std::swap(inp[i], inp[j]);
        }
        else{
            local_high_len++;
        }
    }

    // All-gather on lengths
    // TODO
    int low_len[p], high_len[p];

    // Compute low/high partitions and how many procs to assign for each
    int p_low, p_high;
    int sum_low = 0, sum_high = 0;
    for (int i = 0; i < p; i++) {
        sum_low += low_len[i]; 
        sum_high += high_len[i];
    }
    p_low = std::round(static_cast<float>(sum_low * p) / (sum_high + sum_low));
    // Make sure no problem is left empty
    p_low = std::min(p_low, p - 1);
    p_low = std::max(p_low, 1);
    p_high = p - p_low;

    // Rank belongs to lower half
    int is_lower_half = rank < p_low;
    
    // Split into different communicator
    MPI_Comm my_comm;
    MPI_Comm_split(comm, is_lower_half, rank, &my_comm);

    // Calculate lengths of the new array
    int new_len, new_global_len;
    if(is_lower_half == 1){
        new_global_len = sum_low;
        new_len = sum_low/p_low + (rank < (sum_low % p_low));
    }
    else {
        new_global_len = sum_high;
        new_len = sum_high/p_high + (rank < (sum_high % p_high));
    }

    // Compute displacements, send and receive based on lower half or higher half
    int *sdispls = new int[p], *rdispls = new int[p]; 
    int *scounts = new int[p], *rcounts = new int[p];
    // scount
    for(int i = 0; i < p; i++){
        // Send low elements
        if(i < p_low){
            int low_avail_start = prefix_sum_low[rank];
            int low_avail_end = low_avail_start + local_low_len - 1;
            int low_req_start = prefix_sum_elements[i];
            int low_req_end = low_req_start + num_elements[i];
            scounts[i] = count_overlaps(low_avail_start, low_avail_end, low_req_start, low_req_end);
        }
        // Send high elements
        else {
            int high_avail_start = prefix_sum_high[rank];
            int high_avail_end = high_avail_start + local_high_len - 1;
            int high_req_start = prefix_sum_elements[i] - sum_low;
            int high_req_end = high_req_start + num_elements[i];
            scounts[i] = count_overlaps(high_avail_start, high_avail_end, high_req_start, high_req_end);
        }
    }
    // rcount
    for(int i = 0; i < p; i++){

    }
    
    // All-to-all Communication to split the data into high and low
    MPI_Alltoallv(inp, scounts, sdispls, MPI_INT, inp, rcounts, rdispls, MPI_INT, comm);

    // Compute new seeds
    int seed_low, seed_high;

    if(new_len){
        parallel_qsort(inp, new_len, new_global_len, seed_low, my_comm);
    }
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

