#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <mpi.h>

#define EPS 1e-9
#define MAX_ITER 1000000

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


int main(int argc, char *argv[])
{
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
    int periods[2] = {1, 1};
    int dims[2] = {grid_info.grid_size, grid_info.grid_size};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_info.grid_comm);
    MPI_Cart_coords(grid_info.grid_comm, rank, 2, grid_info.grid_coords);

    std::cout << "Rank " << rank << " coords: " << grid_info.grid_coords[0] << ", " << grid_info.grid_coords[1] << std::endl;
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
    g.global_n = 256;
    int local_ni = 16;
    int local_nj = 16;
    A.resize(local_ni);
    for (int i = 0; i < local_ni; i++)
        A[i].resize(local_nj);
    if (g.grid_coords[1] == 0)
        b.resize(local_ni);
}


void gather_output(char *op_fname, const Vec &x, const GridInfo &g)
{
    // TODO
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
    // Create x_t on all processes
    Vec x_t(A.size());
    // 1. Send x to diagonal
    {
        if (g.grid_coords[0] == 0 && g.grid_coords[1] == 0)
        {
            // Copy
            x_t = x;
            std::cout << g2s(g) << " copied x" << std::endl;
        }
        else if (g.grid_coords[1] == 0)
        {
            // Send
            int recv_rank;
            int recv_coord[2] = {g.grid_coords[0], g.grid_coords[0]};
            MPI_Cart_rank(g.grid_comm, recv_coord, &recv_rank);
            std::cout << g2s(g) << " send to " << recv_rank << std::endl;
            MPI_Send(&x[0], x.size(), MPI_DOUBLE, recv_rank, 0, g.grid_comm);
        }
        else if (g.grid_coords[0] == g.grid_coords[1])
        {
            // Receive
            int send_rank;
            int send_coord[2] = {g.grid_coords[0], 0};
            MPI_Cart_rank(g.grid_comm, send_coord, &send_rank);
            std::cout << g2s(g) << " receive from " << send_rank << std::endl;
            MPI_Recv(&x_t[0], x_t.size(), MPI_DOUBLE, send_rank, 0, g.grid_comm, MPI_STATUS_IGNORE);
        }
    }
    // 2. Broadcast from diagonal
    {
        // Split grid communicator into column communicator
        MPI_Comm col_comm;
        int remain_dims[2] = {0, 1}; // keep column dimension
        MPI_Cart_sub(g.grid_comm, remain_dims, &col_comm);
        // Get rank of diagonal
        // Broadcast 
        int diag_rank;
        int diag_coord = g.grid_coords[0];
        MPI_Cart_rank(col_comm, &diag_coord, &diag_rank);
        // Broadcast 
        MPI_Bcast(&x_t[0], x_t.size(), MPI_DOUBLE, diag_rank, col_comm);
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
    int remain_dims[2] = {1, 0}; // keep row dimension
    MPI_Cart_sub(g.grid_comm, remain_dims, &row_comm);
    // Find rank of column 0
    int root_rank;
    int root_coord = 0;
    MPI_Cart_rank(row_comm, &root_coord, &root_rank);
    // Reduce
    MPI_Reduce(&local_y[0], &y[0], y.size(), MPI_DOUBLE, MPI_SUM, root_rank, row_comm);
}


/*
 * Computes a = a - b
 */
void inplace_vec_sub(Vec &a, const Vec &b)
{
    // TODO
}
