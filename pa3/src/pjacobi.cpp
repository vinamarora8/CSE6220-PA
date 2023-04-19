#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>
#include <mpi.h>

#define EPS 1e-9
#define MAX_ITER 1000000
#define ROOT 0

#ifdef NDEBUG
#   define DBGMSG(en, msg)
#else
#   define DBGMSG(en, msg) {if (en) (std::cout << msg << std::endl);}
#endif

// Types
typedef std::vector<std::vector<double>> Mat;
typedef std::vector<double> Vec;
typedef struct {
    // Original topology
    int rank;
    int world_size;

    // Cartesian topology
    MPI_Comm grid_comm;
    int grid_coords[2]; // [0] = row, [1] = col
    int grid_size; // Problem assumes square grid

    int global_n; // Global size of matrix
} GridInfo;

// Function prototypes
void distribute_inp(Mat &A, Vec &b, GridInfo &g, const char *mat_fname, const char *vec_fname);
void gather_output(char *op_fname, const Vec &x, const GridInfo &g);
void pjacobi_iteration(Vec &x, const Mat &A, const Vec &b, const Vec &d, const GridInfo &g);
double compute_error(const Mat &A, const Vec &x, const Vec &b, const GridInfo &g);
void mat_vec_mult(Vec &y, const Mat &A, const Vec &x, const GridInfo &g, bool ign_diag = false);
void inplace_vec_sub(Vec &a, const Vec &b);
void compute_diagonal(Vec &d, const Mat &A, const GridInfo &g);

// Debug function prototypes
std::string g2s(const GridInfo &g); // Converts grid coords to string


int main(int argc, char *argv[])
{
    bool debug = true;

    char *in_mat_fname = argv[1];
    char *in_vec_fname = argv[2];
    char *op_fname = argv[3];

    int world_size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create grid communicator
    GridInfo grid_info;
    grid_info.rank = rank;
    grid_info.world_size = world_size;
    grid_info.grid_size = sqrt(world_size);
    assert(grid_info.grid_size * grid_info.grid_size == world_size);
    int periods[2] = {0, 0};
    int dims[2] = {grid_info.grid_size, grid_info.grid_size};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_info.grid_comm);
    MPI_Cart_coords(grid_info.grid_comm, rank, 2, grid_info.grid_coords);

    DBGMSG(debug, "Rank " << rank << " coords: " << g2s(grid_info));
    MPI_Barrier(grid_info.grid_comm);

    // Read input matrix and vector
    Mat A;
    Vec b;
    distribute_inp(A, b, grid_info, in_mat_fname, in_vec_fname);
    // // grid_info.global_n, A, b, are now set

    double starttime = MPI_Wtime();

    // Compute answer
    Vec x(b.size());
    Vec d(b.size());
    compute_diagonal(d, A, grid_info);
    double error = compute_error(A, x, b, grid_info);
    int iter = 0;
    while (error > EPS && iter < MAX_ITER)
    {
        pjacobi_iteration(x, A, b, d, grid_info);
        error = compute_error(A, x, b, grid_info);
        iter++;
    }
    // x is now set

    double runtime = (MPI_Wtime() - starttime) * 1000.0;
    if (rank == 0)
        std::cout << "Runtime: " << runtime << " ms" << std::endl;

    // Write output vector
    gather_output(op_fname, x, grid_info);
    
    MPI_Finalize();
}


/*
 * Converts a GridInfo.grid_coords to a string
 */
std::string g2s(const GridInfo &g)
{
    std::stringstream ss;
    ss << "(" << g.grid_coords[0] << ", " << g.grid_coords[1] << ")";
    return ss.str();
}

void distribute_inp(Mat &A, Vec &b, GridInfo &g, const char *mat_fname, const char *vec_fname)
{
    bool debug = false;
    int rank, size, q, n;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    q = sqrt(size);
    if (q * q != size)
    {
        if (rank == 0)
            std::cout << "Number of processes must be a perfect square" << std::endl;
        MPI_Finalize();
        exit(1);
    }

    // Set values for grid info for the given rank of the processor

    if (rank == ROOT)
    {
        // Read the file
        std::ifstream mat_file(mat_fname);
        mat_file >> n;  // Read the size of the matrix
        g.global_n = n; // Set the global size of the matrix
    }
    
    // Broadcast the size of the matrix
    MPI_Bcast(&n, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    // initialize the matrix and vector
    std::vector<int> counts, displs;        // counts and displs for b
    std::vector<int> fill_counts, fill_displs; // counts and displs for A
    std::vector<double> inp, inp_b;         // input buffer for A and b
    int local_ni, local_nj;


    int n_by_q = n / q; // floor(n/q)
    int n_mod_q = n % q;

    // counts, displs, fill_counts, fill_displs are calculated on all the processors since communicating them would be expensive
    counts.resize(q);
    displs.resize(q);

    for (int i = 0; i < q; i++)
    {
        counts[i] = (n_by_q + (i < n_mod_q)); // ceil(n/q) or floor(n/q)
        displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
    }

    // count the number of elements to be filled per processor
    fill_counts.resize(size);
    fill_displs.resize(size);
    
    for (int i = 0; i < size; i++)
    {
        local_ni = counts[i / q];
        local_nj = counts[i % q];
        fill_counts[i] = local_ni * local_nj;
        fill_displs[i] = (i == 0) ? 0 : fill_displs[i - 1] + fill_counts[i - 1];
    }

    if (rank == ROOT)
    {        
        inp.resize(n*n);
        // Load the inp file into matrix inp in a way its continous in memory for a given processor
        double val;
        int count_i, count_j;
        std::ifstream mat_file(mat_fname);
        double tmp;
        mat_file >> tmp; // Need to read again

        for (int i = 0; i < q; i++)
        {
            count_i = counts[i];
            for (int k = 0; k < count_i; k++)
            {
                for (int j = 0; j < q; j++)
                {
                    count_j = counts[j];
                    for (int l = 0; l < count_j; l++)
                    {
                        mat_file >> val;
                        int p = i * q + j;
                        int index = k * count_j + l;
                        inp[fill_displs[p] + index] = val;
                    }
                }
            }
        }

        // Load the vector file into vector b
        std::ifstream vec_file(vec_fname);
        inp_b.resize(n);

        for (int i = 0; i < n; i++)
        {
            double val;
            vec_file >> val;
            inp_b[i] = val;
        }
    }

    local_ni = counts[int(rank / q)];
    local_nj = counts[int(rank % q)];

    A.resize(local_ni);
    for (int i = 0; i < local_ni; i++)
    {
        A[i].resize(local_nj);
    }
    
    if (g.grid_coords[1] == 0)
    {
        b.resize(local_ni);
    }

    // for p processors, scatter the matrix and vector to the processors
    std::vector<double> recvbuff(fill_counts[rank]);
    MPI_Scatterv(
        (void *) &inp[0], (int *) &fill_counts[0], (int *) &fill_displs[0], MPI_DOUBLE,
        (void *) &recvbuff[0], fill_counts[rank], MPI_DOUBLE, ROOT, 
        MPI_COMM_WORLD
    );

    for ( int i = 0; i < local_ni; i++ )
    {
        for ( int j = 0; j < local_nj; j++ )
        {
            A[i][j] = recvbuff[i * local_nj + j];
        }
    }
    // A is set

    // Split communicator for b
    MPI_Comm my_comm;
    MPI_Comm_split(MPI_COMM_WORLD, (rank % q == 0) , rank, &my_comm);

    if (g.grid_coords[1] == 0)
    {
        MPI_Scatterv(
            (void *) &inp_b[0], (int *) &counts[0], (int *) &displs[0], MPI_DOUBLE, 
            (void *) &b[0], counts[int(rank/q)], MPI_DOUBLE, ROOT, 
            my_comm
        );
    }
    // b is set

}

void gather_output(char *op_fname, const Vec &x, const GridInfo &g)
{
    bool debug = false;
    int rank, size, q, n;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (rank == ROOT){
        n = g.global_n;
    }
    MPI_Bcast(&n, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    q = int(sqrt(size));

    // same count and displs can be used for row and column
    int n_by_q = n / q; // floor(n/q)
    int n_mod_q = n % q;
    std::vector<int> counts, displs;
    counts.resize(q);
    displs.resize(q);


    for (int i = 0; i < q; i++)
    {
        counts[i] = (n_by_q + (i < n_mod_q)); // ceil(n/q) or floor(n/q)
        displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
    }
    
    // Gather outputs
    int local_index = int(rank / q);
    Vec comb_array;
    comb_array.resize(n);
    MPI_Comm my_comm;
    MPI_Comm_split(MPI_COMM_WORLD, (rank % q == 0) , rank, &my_comm);

    if (g.grid_coords[1] == 0)
    {

        MPI_Gatherv(
            (void *) &x[0], counts[local_index], MPI_DOUBLE,
            (void *) &comb_array[0], (int *) &counts[0], (int *) &displs[0], 
            MPI_DOUBLE, ROOT, my_comm
        );
    }

    if (rank == ROOT)
    {
        std::ofstream opfile;
        opfile.open(op_fname);
        opfile << std::fixed << std::setprecision(16);
        for (int i = 0; i < n; i++)
        {
            opfile << comb_array[i] << " ";
        }
        opfile << std::endl;
        opfile.close();
    }

}

/*
 * Computes d = Diagonal vec of A
 */
void compute_diagonal(Vec &d, const Mat &A, const GridInfo &g)
{
    // compute local_d
    Vec local_d(A.size());
    for(int i = 0; i < local_d.size(); i++)
    {
        local_d[i] = A[i][i];
    }
    // Send
    if(g.grid_coords[0] == g.grid_coords[1]){
        int recv_rank;
        int recv_coord[2] = {g.grid_coords[0], 0};
        MPI_Cart_rank(g.grid_comm, recv_coord, &recv_rank);
        MPI_Send(&local_d[0], local_d.size(), MPI_DOUBLE, recv_rank, 0, g.grid_comm);
    }
    // Receive
    if(g.grid_coords[1] == 0){
        int send_rank;
        int send_coord[2] = {g.grid_coords[0], g.grid_coords[0]};
        MPI_Cart_rank(g.grid_comm, send_coord, &send_rank);
        MPI_Recv(&d[0], local_d.size(), MPI_DOUBLE, send_rank, 0, g.grid_comm, MPI_STATUS_IGNORE);
    }
}

/*
 * Computes x = D^-1 (b - Rx)
 */
void pjacobi_iteration(Vec &x, const Mat &A, const Vec &b, const Vec &d, const GridInfo &g)
{
    Vec y(A[0].size());
    // Compute Rx
    mat_vec_mult(y, A, x, g, true);
    // Compute b - Rx
    inplace_vec_sub(y, b);
    // Compute x = D_inv(b - Rx) for column 0
    if(g.grid_coords[1] == 0)
    {
        for(int i = 0; i < A.size(); i++)
        {
            y[i] = y[i] / d[i];
        }
    }
    x = y;
}


/*
 * Computes L2 error between Ax and b
 */
double compute_error(const Mat &A, const Vec &x, const Vec &b, const GridInfo &g)
{
    double err = 0.0;
    double total_err = 0.0;
    Vec y;
    
    mat_vec_mult(y, A, x, g, false); // need to consider diagonal elements
    
    // TODO changes pertaining to b being a column vector
    for(int i = 0; i < y.size(); i++){
        err += ((y[i] - b[i]) * (y[i] - b[i]));
    }
    
    // MPI_Reduce(&err, &total_err, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
    // all reduce
    MPI_Allreduce(&err, &total_err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    total_err = sqrt(total_err);

    return total_err;
}


/*
 * Computes y = Ax
 * If ign_diag is true, then the diagonal elements of A are ignored (taken to be 0)
 */
void mat_vec_mult(Vec &y, const Mat &A, const Vec &x, const GridInfo &g, bool ign_diag)
{
    bool debug = false;

    // Create x_t on all processes.
    Vec x_t(A[0].size());
    // Two steps:
    // 1. Send x to diagonal
    {
        if (g.grid_coords[0] == 0 && g.grid_coords[1] == 0)
        {
            // Copy
            // Separate so that Send and Recv don't overlap for the same process
            x_t = x;
            DBGMSG(debug, g2s(g) << " copied x");
        }
        else if (g.grid_coords[1] == 0)
        {
            // Send
            int recv_rank;
            int recv_coord[2] = {g.grid_coords[0], g.grid_coords[0]};
            MPI_Cart_rank(g.grid_comm, recv_coord, &recv_rank);
            MPI_Send(&x[0], x.size(), MPI_DOUBLE, recv_rank, 0, g.grid_comm);
            DBGMSG(debug, g2s(g) << " send to " << recv_rank);
        }
        else if (g.grid_coords[0] == g.grid_coords[1])
        {
            // Receive
            int send_rank;
            int send_coord[2] = {g.grid_coords[0], 0};
            MPI_Cart_rank(g.grid_comm, send_coord, &send_rank);
            MPI_Recv(&x_t[0], x_t.size(), MPI_DOUBLE, send_rank, 0, g.grid_comm, MPI_STATUS_IGNORE);
            DBGMSG(debug, g2s(g) << " receive from " << send_rank);
        }
    }
    // 2. Broadcast from diagonal
    {
        // Split grid communicator into column communicator
        MPI_Comm col_comm;
        int remain_dims[2] = {1, 0}; // drop column dimension
        MPI_Cart_sub(g.grid_comm, remain_dims, &col_comm);
        // Find rank of diagonal
        int diag_rank;
        int diag_coord[1] = {g.grid_coords[1]};
        MPI_Cart_rank(col_comm, diag_coord, &diag_rank);
        // Broadcast 
        MPI_Bcast(&x_t[0], x_t.size(), MPI_DOUBLE, diag_rank, col_comm);
        DBGMSG(debug, g2s(g) << " got x_t from " << diag_rank);
    }

    // Local dot product
    Vec local_y(A.size());
    for (int i = 0; i < A.size(); i++)
    {
        double sum = 0.0;
        for (int j = 0; j < A[i].size(); j++)
        {
            if (ign_diag && g.grid_coords[0] == g.grid_coords[1] && i == j)
                continue;
            sum += A[i][j] * x_t[j];
        }
        local_y[i] = sum;
    }

    // Reduce local_y to y
    // Split grid communicator into row communicator
    MPI_Comm row_comm;
    int remain_dims[2] = {0, 1}; // drop row dimension 
    MPI_Cart_sub(g.grid_comm, remain_dims, &row_comm);
    // Find rank of column 0
    int root_rank;
    int root_coord[1] = {0};
    MPI_Cart_rank(row_comm, root_coord, &root_rank);
    // Reduce
    MPI_Reduce(&local_y[0], &y[0], local_y.size(), MPI_DOUBLE, MPI_SUM, root_rank, row_comm);
    DBGMSG(debug, g2s(g) << " reduced to " << root_rank);
}


/*
 * Computes a = b - a
 */
void inplace_vec_sub(Vec &a, const Vec &b)
{
    // TODO
    for (int i = 0; i < a.size(); i++)
    {
        a[i] = b[i] - a[i];
    }
}
