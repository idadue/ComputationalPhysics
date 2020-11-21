#include "ising.h"
#include "mpi.h"
#include <time.h>

//TODO Add timing

int main(int argc, char *argv[])
{
    clock_t start, finish;

    int rank, processes;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double init_temp = 2.0;
    double final_temp = 2.3;
    double dT = 0.05;
    int cycles = 1000000;
    int L = 60;

    int no_intervals = cycles / processes;
    int myloop_begin = rank * no_intervals + 1;
    int myloop_end = (rank + 1) * no_intervals;

    if ((rank == processes - 1) && (myloop_end < cycles))
        myloop_end = cycles;

    // broadcast to all nodes common variables
    MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&init_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    Ising is(L, "test_mpi2");
    is.initializeSystem(init_temp);
    is.print();
    std::cout << "------------\n";
    start = clock();
    for (double temp = init_temp; temp <= final_temp; temp += dT)
    {
        is.updateSystem(temp);
        is.MonteCarloCycler(cycles);
    }
    finish = clock();
    double elapsed_time = double(finish - start) / CLOCKS_PER_SEC;
    is.print();
    printf("t = %.3f \bs\n", elapsed_time);

    MPI_Finalize();

    return 0;
}