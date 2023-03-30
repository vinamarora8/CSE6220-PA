#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <string>

int main(int argc, char* argv[])
{
    int rank, size;
    int n;

    std::cout << std::fixed << std::setprecision(12);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Process cl args
    if (rank == 0)
    {
        if (argc < 2)
        {
            std::cerr << "Missing argument: n" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 255);
        }
        else
        {
            n = std::stoi(argv[1]);
        }
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double starttime = MPI_Wtime();

    // Compute work boundaries
    int start, end;
    start = (n * rank / size) + 1;
    end = n * (rank + 1) / size;

    // Local sums
    double local_sum = 0.0, h;
    for (int i = start; i <= end; i++)
    {
        h = ((double) i - 0.5) / (double) n;
        h *= h;
        h = 1.0 / (1.0 + h);
        local_sum += h;
    }

    // Reduction
    double sum;
    MPI_Reduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Output
    double runtime = MPI_Wtime() - starttime;
    if (rank == 0)
    {
        sum *= 4.0 / (double) n;
        std::cout << sum << ", " << runtime << std::endl;
    }

    MPI_Finalize();
    return 0;
}
