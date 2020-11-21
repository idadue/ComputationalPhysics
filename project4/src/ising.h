#include <iostream>
#include "random.h"
#include <fstream>
#include <iomanip>

class Ising
{
    /*
    Class for simulating the two-dimensional Ising model. We define 
    the coupling constant as J = 1 and the Boltzmann constant as k_B = 1.
    Then the temperature has the dimesion [Energy]. Simulation is done using 
    Monte Carlo cylces via the Metropolis algortihm with periodic 
    boundary conditions.    
    */

public:
    Ising(int L, std::string filename);
    ~Ising();
    void initializeSystem(double temperature);
    void updateSystem(double temperature);
    void Metropolis();
    void MonteCarloCycler(int cycles);
    void computeValuesOfInterest(int cycle);

    void writeToFile(std::string filename);

    void print()
    {
        for (int row = 0; row < L; row++)
        {
            for (int col = 0; col < L; col++)
            {
                printf("%3d ", spin_matrix[row * L + col]);
            }
            std::cout << "\n";
        }
    };

private:
    int L;                 //Square lattice dimension
    int N;                 //Number of spins
    int cycles;            //Number of Monte Carlo cycles to run
    double w[17];          //energy difference
    double expectation[5]; //expectation values

    double E;
    double E_mean = 0;

    double M;
    double M_mean = 0;
    double M_mean_abs = 0;

    double temperature;
    double specific_heat = 0;
    double susceptibility = 0;

    Random r; //Object of Random class
    std::string filename;
    int configurations = pow(2, N);
    int *spin_matrix;

    int periodic(int i, int add);
    int getDeltaE(int row, int col);
    void energyDiff(double temperature);
};