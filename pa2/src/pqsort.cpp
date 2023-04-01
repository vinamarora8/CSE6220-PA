#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>
#include <mpi.h>

#define ROOT 0
//#define DEBUG_BLOCK_DIST
//#define DEBUG_ALLTOALL
//#define DEBUG_PIVOT

void serial_sort(int *inp, int low, int high);
void parallel_qsort(int *&inp, int &len, int global_len, int seed, MPI_Comm comm);
int distribute_input(const char *fname, int *&local_inp, int *global_len, MPI_Comm comm);
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
    int local_len, *local_inp, global_len;
    local_len = distribute_input(in_fname, local_inp, &global_len, MPI_COMM_WORLD);

    // Timing start
    double starttime = MPI_Wtime();

    parallel_qsort(local_inp, local_len, global_len, 0, MPI_COMM_WORLD);

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

// Computes the number of overlaps between two ranges.
int count_overlaps(int start1, int end1, int start2, int end2) {
    int max_start = std::max(start1, start2);
    int min_end = std::min(end1, end2);
    int overlap_count = std::max(0, min_end - max_start + 1);
    return overlap_count;
}
// Computes the exclusive prefix sum of an array.
void exclusive_prefix_sum(const int* in, int* out, int size) {
    int sum = 0;
    for (int i = 0; i < size; i++) {
        out[i] = sum;
        sum += in[i];
    }
}

/**
 * Compute lowest global_index in rank given world size p and array length
*/
int global_index_low(int global_len, int p, int rank)
{
    int ans = rank * (global_len / p);

    if (rank < (global_len % p))
        ans += rank;
    else
        ans += global_len % p;
    
    return ans;
}
std::string arrayToString(int* arr, int size) {
  std::stringstream ss;
  ss << "[";
  for (int i = 0; i < size; i++) {
    ss << arr[i];
    if (i < size - 1) {
      ss << ", ";
    }
  }
  ss << "]";
  return ss.str();
}

void parallel_qsort(int *&inp, int &len, int global_len, int seed, MPI_Comm comm)
{

    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    if (p == 1)
    {
        serial_sort(inp, 0, len-1);
        //std::printf("Rank:%d, inp: %s\n", rank, arrayToString(inp, len).c_str());
        return;
    }

    // Choose pivot (random, with same seed) and broadcast
    std::srand(seed);

    int pivot_index = std::rand() % global_len;
    int pivot = -1;

    // find which processor has the pivot, broadcast to all
    int pivot_rank;
    for (pivot_rank = 0; pivot_rank < p; pivot_rank++)
    {
        if (global_index_low(global_len, p, pivot_rank) <= pivot_index 
            && global_index_low(global_len, p, pivot_rank+1) > pivot_index)
            break;
    }
    if (rank == pivot_rank)
        pivot = inp[pivot_index - global_index_low(global_len, p, rank)];
    MPI_Bcast(&pivot, 1, MPI_INT, pivot_rank, comm);

#ifdef DEBUG_PIVOT
    if (rank == 0)
    {
        std::cout << "Rank:0, Pivot_idx: " << pivot_index << std::endl;
        std::cout << "Rank:0, Pivot_rank: " << pivot_rank << std::endl;
        std::cout << "Rank:0, Pivot: " << pivot << std::endl;
    }
#endif


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
    int* low_len = new int[p];
    int* high_len = new int[p];

    MPI_Allgather(&local_low_len, 1, MPI_INT, low_len, 1, MPI_INT, comm);
    MPI_Allgather(&local_high_len, 1, MPI_INT, high_len, 1, MPI_INT, comm);

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
        new_len = sum_high/p_high + ((rank-p_low) < (sum_high % p_high));
    }

    // Num elements
    int *num_elements = new int[p];
    for (int i = 0; i < p; i++) {
        if(i < p_low){
            num_elements[i] = sum_low/p_low + (i < (sum_low % p_low));
        }
        else{
            num_elements[i] = sum_high/p_high + ((i-p_low) < (sum_high % p_high));
        }
    }
    // Prefix sum calculations required
    int *prefix_sum_elements = new int[p];
    int *prefix_sum_low = new int[p];
    int *prefix_sum_high = new int[p];
    exclusive_prefix_sum(num_elements, prefix_sum_elements, p);
    exclusive_prefix_sum(low_len, prefix_sum_low, p);
    exclusive_prefix_sum(high_len, prefix_sum_high, p);

    // Compute displacements, send and receive based on lower half or higher half
    int *sdispls = new int[p], *rdispls = new int[p]; 
    int *scounts = new int[p], *rcounts = new int[p];
    // scount
    for(int i = 0; i < p; i++){
        scounts[i] = 0;
        // Send low elements
        if(i < p_low){
            int low_avail_start = prefix_sum_low[rank];
            int low_avail_end = low_avail_start + local_low_len - 1;
            int low_req_start = prefix_sum_elements[i];
            int low_req_end = low_req_start + num_elements[i] - 1;
            if(low_req_start >= 0 && low_req_end >= 0 && low_avail_start >= 0 && low_avail_end >= 0){
                scounts[i] = count_overlaps(low_avail_start, low_avail_end, low_req_start, low_req_end);
            }
        }
        // Send high elements
        else {
            int high_avail_start = prefix_sum_high[rank];
            int high_avail_end = high_avail_start + local_high_len - 1;
            int high_req_start = prefix_sum_elements[i] - sum_low;
            int high_req_end = high_req_start + num_elements[i] - 1;
            if(high_req_start >= 0 && high_req_end >= 0 && high_avail_start >= 0 && high_avail_end >= 0){
                scounts[i] = count_overlaps(high_avail_start, high_avail_end, high_req_start, high_req_end);
            }
        }
    }
    // rcount
    for(int i = 0; i < p; i++){
        rcounts[i] = 0;
        // Receive low elements
        if(rank < p_low){
            int low_req_start = prefix_sum_elements[rank];
            int low_req_end = low_req_start + num_elements[rank] - 1;
            int low_avail_start = prefix_sum_low[i];
            int low_avail_end = low_avail_start + low_len[i] - 1;
            if(low_req_start >= 0 && low_req_end >= 0 && low_avail_start >= 0 && low_avail_end >= 0){
                rcounts[i] = count_overlaps(low_avail_start, low_avail_end, low_req_start, low_req_end);
            }
        }
        else {
            int high_req_start = prefix_sum_elements[rank] - sum_low;
            int high_req_end = high_req_start + num_elements[rank] - 1;
            int high_avail_start = prefix_sum_high[i];
            int high_avail_end = high_avail_start + high_len[i] - 1;
            if(high_req_start >= 0 && high_req_end >= 0 && high_avail_start >= 0 && high_avail_end >= 0){
                rcounts[i] = count_overlaps(high_avail_start, high_avail_end, high_req_start, high_req_end);
            }
        }
    }
    // sdispls
    int sum_displ = 0;
    for(int i = 0; i < p; i++){
        sdispls[i] = sum_displ;
        sum_displ += scounts[i]; 
    }    
    // rdispls
    sum_displ = 0;
    for(int i = 0; i < p; i++){
        rdispls[i] = sum_displ;
        sum_displ += rcounts[i]; 
    }

#ifdef DEBUG_ALLTOALL
    std::printf("Rank is %d, pivot is %d, p_low: %d, p_high: %d\n", rank, pivot, p_low, p_high);
    std::printf("Rank:%d, inp: %s\n", rank, arrayToString(inp, len).c_str());
    std::printf("Rank:%d, low_len: %s\n", rank, arrayToString(low_len, p).c_str());
    std::printf("Rank:%d, high_len: %s\n", rank, arrayToString(high_len, p).c_str());
    std::printf("Rank:%d, prefix_sum_low: %s\n", rank, arrayToString(prefix_sum_low, len).c_str());
    std::printf("Rank:%d, prefix_sum_high: %s\n", rank, arrayToString(prefix_sum_high, len).c_str());
    std::printf("Rank:%d, num_elements: %s\n", rank, arrayToString(num_elements, p).c_str());
    std::printf("Rank:%d, scounts: %s\n", rank, arrayToString(scounts, p).c_str());
    std::printf("Rank:%d, rcounts: %s\n", rank, arrayToString(rcounts, p).c_str());
    std::printf("Rank:%d, sdispls: %s\n", rank, arrayToString(sdispls, p).c_str());
    std::printf("Rank:%d, rdispls: %s\n", rank, arrayToString(rdispls, p).c_str());
#endif

    // All-to-all Communication to split the data into high and low
    int *new_inp = new int[new_len];
    MPI_Alltoallv(inp, scounts, sdispls, MPI_INT, new_inp, rcounts, rdispls, MPI_INT, comm);

#ifdef DEBUG_ALLTOALL
    std::printf("Rank:%d, after_inp: %s\n", rank, arrayToString(new_inp, new_len).c_str());
#endif

    free(inp);
    delete[] rcounts;
    delete[] rdispls;
    delete[] scounts;
    delete[] sdispls;

    inp = new_inp;
    len = new_len;

    if(new_len){
        parallel_qsort(inp, len, new_global_len, std::rand(), my_comm);
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
int distribute_input(const char *fname, int *&local_inp, int *global_len, MPI_Comm comm)
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

    // Broadcast global length
    *global_len = len;
    MPI_Bcast(global_len, 1, MPI_INT, ROOT, comm);

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
            //std::cout << comb_array[i] << " ";
            fstream << comb_array[i] << " ";
        }
        //std::cout << std::endl;
        fstream << std::endl;
    }

    free(comb_array);
}

