#include <iostream>
#include "random.h"
#include <fstream>
#include <iomanip>
#include <vector>

class Ising
{
    /*
    Class for simulating the two-dimensional Ising model. We define 
    the coupling constant as J = 1 and the Boltzmann constant as k_B = 1.
    Simulation is done using Monte Carlo cylces via 
    the Metropolis algortihm with periodic boundary conditions.    
    */

public:
    Ising(int L, std::string filename);
    ~Ising();
    void initializeSystem(double temperature);
    void updateSystem(double temperature);
    void Metropolis();
    void MonteCarloCycler(int cycles, std::vector<double> &expectation);
    void computeValuesOfInterest(int cycle, std::vector<double> &expectation);

    void writeToFile(int cycle, std::vector<double> &expectation);

private:
    int L;        //Square lattice dimension
    int N;        //Number of spins
    int cycles;   //Number of Monte Carlo cycles to run
    double w[17]; //energy difference

    double E;          //energy
    double E_mean = 0; //expectation value of norm of energy

    double M;              //magnetization
    double M_mean = 0;     //expectation value of norm of magnetization
    double M_mean_abs = 0; //absolute value of norm of magnetization

    double temperature;
    double specific_heat = 0;
    double susceptibility = 0;

    Random r; //Object of Random class
    std::string filename;

    int *spin_matrix;

    int periodic(int i, int add);
    int getDeltaE(int row, int col);
    void energyDiff(double temperature);
};