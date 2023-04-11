#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>

#define EPS 1e-9
#define MAX_ITER 1000000

// Global vars that all functions need
double **local_mat, *local_vec, *local_res;
int local_ni, local_nj, global_n;
MPI_Comm grid_comm;

// Problem specific functions
void distribute_inp(char *mat_fname, char *vec_fname);
void gather_output(char *op_fname);
void pjacobi_iteration(double **mat, double *vec);
double compute_error(double **mat, double *vec, double *res);
void mat_vec_mult(double **mat, double *vec, double *res); // size is same as global vars

int main(int argc, char *argv[])
{
    char *in_mat_fname = argv[1];
    char *in_vec_fname = argv[2];
    char *op_fname = argv[3];

    int p, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create grid communicator
    // TODO
    // grid_comm is now set

    // Read input matrix and vector
    distribute_inp(in_mat_fname, in_vec_fname);
    // global_n, local_ni, local_nj, local_mat, local_vec are now set

    double starttime = MPI_Wtime();

    // Compute answer
    local_res = new double[local_ni];
    double error = 1.0;
    int iter = 0;
    while (error > EPS && iter < MAX_ITER)
    {
        pjacobi_iteration(local_mat, local_vec);
        error = compute_error(local_mat, local_vec, local_res);
        iter++;
    }
    // local_res is now set

    double runtime = (MPI_Wtime() - starttime) * 1000.0;
    if (rank == 0)
        std::cout << "Runtime: " << runtime << " ms" << std::endl;

    // Write output vector
    gather_output(op_fname);
    
    // Free memory
    delete[] local_res;
    delete[] local_vec;
    for (int i = 0; i < local_ni; i++)
        delete[] local_mat[i];

    MPI_Finalize();
}


void distribute_inp(char *mat_fname, char *vec_fname)
{
    // TODO: Have to set these global vars
    global_n = 0;
    local_ni = 0;
    local_nj = 0;
    for (int i = 0; i < local_ni; i++)
        local_mat[i] = new double[local_nj];
    local_vec = new double[local_ni];
}


void gather_output(char *op_fname)
{
    // TODO: Take output from local_res
}


void pjacobi_iteration(double **mat, double *vec)
{
    // TODO: Have to set local_res
}


double compute_error(double **mat, double *vec, double *res)
{
    // TODO
    double err = 0.0;

    return err;
}


void mat_vec_mult(double **mat, double *vec, double *res)
{
    // TODO
}
