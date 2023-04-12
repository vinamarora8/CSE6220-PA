#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <mpi.h>

#define EPS 1e-9
#define MAX_ITER 1000000

// Types
typedef std::vector<std::vector<double>> Mat;
typedef std::vector<double> Vec;
typedef struct {
    MPI_Comm comm;
    int coords[2]; // [0] = row, [1] = col
    int size; // Problem assumes square grid
    int global_n; // Global size of matrix
} GridInfo;

// Function prototypes
void distribute_inp(Mat &A, Vec &b, GridInfo &g, const char *mat_fname, const char *vec_fname);
void gather_output(char *op_fname, const Vec &x, const GridInfo &g);
void pjacobi_iteration(Vec &x, const Mat &A, const Vec &b, const GridInfo &g);
double compute_error(const Mat &A, const Vec &x, const Vec &b, const GridInfo &g);
void mat_vec_mult(Vec &y, const Mat &A, const Vec &x, const GridInfo &g);
void local_vec_sub(Vec &y, const Vec &x1, const Vec &x2);


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
    // TODO
    GridInfo grid_info;
    grid_info.size = sqrt(world_size);
    // grid_info.comm and grid_info.size are now set

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


void distribute_inp(Mat &A, Vec &b, GridInfo &g, const char *mat_fname, const char *vec_fname)
{
    // TODO: Have to set these and fill the matrix
    int local_ni = 0;
    int local_nj = 0;
    g.global_n = 0;
    A.resize(local_ni);
    for (int i = 0; i < local_ni; i++)
        A[i].resize(local_nj);
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
 */
void mat_vec_mult(Vec &y, const Mat &A, const Vec &x, const GridInfo &g)
{
    // TODO
}


/*
 * Computes y = x1 - x2
 */
void local_vec_sub(Vec &y, const Vec &x1, const Vec &x2, const GridInfo &g)
{
    // TODO
}
