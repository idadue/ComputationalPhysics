#include "ising.h"
#include "mpi.h"
#include <time.h>

int main(int argc, char *argv[])
{
    clock_t start, finish;
    std::string filename;
    double init_temp, final_temp, dT;
    int cycles, L;

    int rank, processes;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0 && argc > 1)
    {
        filename = argv[1];
        L = atoi(argv[2]);
        cycles = atoi(argv[3]);
        init_temp = atof(argv[4]);
        final_temp = atof(argv[5]);
        dT = atof(argv[6]);
    }

    // broadcast to all nodes common variables
    MPI_Bcast(&cycles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&init_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    Ising is(L, filename);
    is.initializeSystem(init_temp);
    start = clock();
    double TimeStart, TimeEnd, TotalTime;
    TimeStart = MPI_Wtime();
    for (double temp = init_temp; temp <= final_temp; temp += dT)
    {
        std::vector<double> localExpectation = {0, 0, 0, 0, 0};
        std::vector<double> globalExpectation = {0, 0, 0, 0, 0};

        is.updateSystem(temp);
        is.MonteCarloCycler(cycles, localExpectation);

        /*for (int i = 0; i < 5; i++)
        {
            MPI_Reduce(&localExpectation[i], &globalExpectation[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }*/
        if (rank == 0)
        {
            //is.writeToFile(cycles, globalExpectation);
            is.writeToFile(cycles, localExpectation);
        }
    }
    finish = clock();
    TimeEnd = MPI_Wtime();
    TotalTime = TimeEnd - TimeStart;
    if (rank == 0)
    {
        std::cout << "Time = " << TotalTime << " on number of processors: " << processes << std::endl;
    }
    double elapsed_time = double(finish - start) / CLOCKS_PER_SEC;
    printf("t = %.3f \bs\n", elapsed_time);

    MPI_Finalize();

    return 0;
}