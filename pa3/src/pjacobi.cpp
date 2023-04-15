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
void pjacobi_iteration(Vec &x, const Mat &A, const Vec &b, const GridInfo &g);
double compute_error(const Mat &A, const Vec &x, const Vec &b, const GridInfo &g);
void mat_vec_mult(Vec &y, const Mat &A, const Vec &x, const GridInfo &g, bool ign_diag = false);
void inplace_vec_sub(Vec &a, const Vec &b);

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
    // grid_info.global_n, A, b, are now set

    double starttime = MPI_Wtime();

    // Compute answer
    Vec x(b.size());
    double error = 1.0;
    int iter = 0;
    while (error > EPS && iter < MAX_ITER)
    {
        pjacobi_iteration(x, A, b, grid_info);
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
    // TODO: Have to set these and fill the matrix
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
    g.grid_coords[0] = rank / q;
    g.grid_coords[1] = rank % q;

    if (rank == ROOT)
    {
        // Read the file
        std::ifstream mat_file(mat_fname);
        mat_file >> n;  // Read the size of the matrix
        g.global_n = n; // Set the global size of the matrix
    }

    MPI_Bcast(n, 1, MPI_INT, ROOT, comm);

    // initialize the matrix and vector
    std:: : vector<int> count, displs;        // counts and displs for b
    std::vector<int> fill_count, fill_displs; // counts and displs for A
    std::vector<double> inp;
    std::vector<double> inp_b;

    // calculate the counts and displacements for the scatterv function
    int local_ni, local_nj;

    // same count and displs can be used for row and column
    int n_by_q = n / q; // floor(n/q)
    int n_mod_q = n % q;

    count.resize(q);
    displs.resize(q);


    for (int i = 0; i < q; i++)
    {
        count[i] = (n_by_q + (i < n_mod_q)); // ceil(n/q) or floor(n/q)
        displs[i] = (i == 0) ? 0 : displs[i - 1] + count[i - 1];
    }

    if (rank == ROOT)
    {
        // count the number of elements to be filled per processor
        fill_count.resize(size);
        fill_displs.resize(size);
        inp.resize(size);

        for (i = 0; i < size; i++)
        {
            local_ni = count[i / q];
            local_nj = count[i % q];
            fill_count[i] = local_ni * local_nj;
            fill_displs[i] = (i == 0) ? 0 : fill_displs[i - 1] + fill_count[i - 1];
            inp[i].resize(fill_count[i]);
        }

        // Load the inp file into matrix inp in a way its continous in memory for a given processor
        double val;
        for (int i = 0; i < q; i++)
        {
            count_i = count[i];
            for (int k = 0; k < count_i; k++)
            {
                for (int j = 0; j < q; j++)
                {
                    count_j = count[j];
                    for (int l = 0; l < count_j; l++)
                    {
                        mat_file >> val;
                        int p = i * q + j;
                        int index = k * count_j + l;
                        inp[p][index] = val;
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
            inp_b[i].push_back(val);
        }
    }

    // set the local_ni and local_nj for the processor
    MPI_Scatter(
        (int *)counts, 1, MPI_INT,
        &local_count, 1, MPI_INT,
        ROOT, comm);

    local_ni = count[rank / q];
    local_ni = count[rank % q];
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
    MPI_Scatterv(
        inp, fill_count, fill_displs, MPI_DOUBLE,
        A, fill_count[rank], MPI_DOUBLE, ROOT, 
        MPI_COMM_WORLD
    );

    MPI_Scatterv(
        inp_b, count, displs, MPI_DOUBLE, 
        b, count[rank], MPI_DOUBLE, ROOT, 
        MPI_COMM_WORLD
    );
}

void gather_output(char *op_fname, const Vec &x, const GridInfo &g)
{
    // TODO
    
    int rank, size, q, n;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    n = g.global_n;

    // same count and displs can be used for row and column
    int n_by_q = n / q; // floor(n/q)
    int n_mod_q = n % q;

    count.resize(q);
    displs.resize(q);

    for (int i = 0; i < q; i++)
    {
        count[i] = (n_by_q + (i < n_mod_q)); // ceil(n/q) or floor(n/q)
        displs[i] = (i == 0) ? 0 : displs[i - 1] + count[i - 1];
    }

    // Gather outputs
    std::vector<double> *comb_array;
    MPI_Gather(
        x, count[rank], MPI_DOUBLE,
        comb_array, counts, displs, MPI_DOUBLE,
        ROOT, MPI_COMM_WORLD
    );

    if (rank == ROOT)
    {
        for (int i = 0; i < total_len; i++)
        {
            op_fname << comb_array[i] << " ";
        }
        op_fname << std::endl;
    }

    free(comb_array);


}


/*
 * Computes x = D^-1 (b - Rx)
 */
void pjacobi_iteration(Vec &x, const Mat &A, const Vec &b, const GridInfo &g)
{
    // TODO
}


/*
 * Computes L2 error between Ax and b
 */
double compute_error(const Mat &A, const Vec &x, const Vec &b, const GridInfo &g)
{
    // TODO
    double err = 0.0;

    return err;
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
        int diag_coord[1] = {g.grid_coords[0]};
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
            if (ign_diag && i == j)
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
 * Computes a = a - b
 */
void inplace_vec_sub(Vec &a, const Vec &b)
{
    // TODO

}
