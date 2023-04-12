#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <mpi.h>

#define EPS 1e-9
#define MAX_ITER 1000000

typedef std::vector<std::vector<double>> MatD;
typedef std::vector<double> VecD;

// Global vars that all functions need
MatD local_mat;
VecD local_vec, local_res;
int global_n;
MPI_Comm grid_comm;
int world_size, rank;
int grid_rank[2], grid_size[2]; // [0] -> rows, [1] -> cols

// Problem specific functions
void distribute_inp(char *mat_fname, char *vec_fname);
void gather_output(char *op_fname);
void pjacobi_iteration(MatD &mat, VecD &vec, VecD &res);
double compute_error(MatD &mat, VecD &vec, VecD &res);
void mat_vec_mult(MatD &mat, VecD &vec, VecD &res);
void local_vec_sub(VecD &vec1, VecD &vec2, VecD &res, int n);


int main(int argc, char *argv[])
{
    char *in_mat_fname = argv[1];
    char *in_vec_fname = argv[2];
    char *op_fname = argv[3];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create grid communicator
    // TODO
    // grid_comm, grid_size, grid_rank is now set

    // Read input matrix and vector
    distribute_inp(in_mat_fname, in_vec_fname);
    // global_n, local_n, local_mat, local_vec are now set

    double starttime = MPI_Wtime();

    // Compute answer
    local_res.resize(local_vec.size());
    double error = 1.0;
    int iter = 0;
    while (error > EPS && iter < MAX_ITER)
    {
        pjacobi_iteration(local_mat, local_vec, local_res);
        error = compute_error(local_mat, local_vec, local_res);
        iter++;
    }
    // local_res is now set

    double runtime = (MPI_Wtime() - starttime) * 1000.0;
    if (rank == 0)
        std::cout << "Runtime: " << runtime << " ms" << std::endl;

    // Write output vector
    gather_output(op_fname);
    
    MPI_Finalize();
}


void distribute_inp(char *mat_fname, char *vec_fname)
{
    // TODO: Have to set these
    int local_ni = 0;
    int local_nj = 0;
    global_n = 0;
    local_mat.resize(local_ni);
    for (int i = 0; i < local_ni; i++)
        local_mat[i].resize(local_nj);
    local_vec.resize(local_ni);
}


void gather_output(char *op_fname)
{
    // TODO: Take output from local_res
}


void pjacobi_iteration(MatD &mat, VecD &vec, VecD &res)
{
    // TODO: Have to set res
}


double compute_error(MatD &mat, VecD &vec, VecD &res)
{
    // TODO
    double err = 0.0;

    return err;
}


void mat_vec_mult(MatD &mat, VecD &vec, VecD &res)
{
    // TODO
}


void local_vec_sub(VecD &vec1, VecD &vec2, VecD &res, int n)
{
    // TODO
}
