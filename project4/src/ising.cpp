#include "ising.h"

Ising::Ising(int L, std::string filename) : L(L), N(L * L), filename(filename)
{
    spin_matrix = new int[N];
    /*
    std::ofstream ofile;

    ofile.open(filename + ".txt", std::ios::out);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << std::setw(15) << "T:";
    ofile << std::setw(15) << "E_mean:";
    ofile << std::setw(15) << "M_mean:";
    ofile << std::setw(15) << "SpH:";
    ofile << std::setw(15) << "S:";
    ofile << std::setw(15) << "M_mean_abs:\n";
    ofile.close();*/

    //initializeSystem();
    //std::cout << "here is fine 1 \n";
}

Ising::~Ising()
{
    delete[] spin_matrix;
}

void Ising::updateSystem(double temperature)
{
    M = 0.0;
    E = 0.0;
    this->temperature = temperature;

    double seed = r.generateMersenneTwisterNumber(0, 10000);
    r.setSeed(seed);

    for (int row = 0; row < L; row++)
    {
        for (int column = 0; column < L; column++)
        {
            if (temperature < 1.5)
            {
                spin_matrix[row * L + column] = 1;
            }

            M += (double)spin_matrix[row * L + column];
        }
    }

    for (int row = 0; row < L; row++)
    {
        for (int column = 0; column < L; column++)
        {
            int p = periodic(row, -1);
            int p2 = periodic(column, -1);
            E -= spin_matrix[row * L + column] * (spin_matrix[p * L + column] + spin_matrix[row * L + p2]);
        }
    }

    energyDiff(temperature);
    for (int i = 0; i < 5; i++)
    {
        expectation[i] = 0;
    }
}

void Ising::initializeSystem(double temperature)
{
    double seed = r.generateMersenneTwisterNumber(0, 10000);
    r.setSeed(seed);
    this->temperature = temperature;

    for (int row = 0; row < L; row++)
    {
        for (int column = 0; column < L; column++)
        {
            spin_matrix[row * L + column] = r.spin();
            //M += (double)spin_matrix[row * L + column];
        }
    }
    /*
    for (int row = 0; row < L; row++)
    {
        for (int column = 0; column < L; column++)
        {
            int p = periodic(row, -1);
            int p2 = periodic(column, -1);
            E -= spin_matrix[row * L + column] * (spin_matrix[p * L + column] + spin_matrix[row * L + p2]);
        }
    }

    energyDiff(temperature);
    for (int i = 0; i < 5; i++)
    {
        expectation[i] = 0;
    }*/
}

int Ising::periodic(int i, int add)
{
    return (i + L + add) % (L);
}

void Ising::energyDiff(double temperature)
{
    for (int i = -8; i <= 8; i++)
    {
        w[i + 8] = 0;
    }
    for (int i = -8; i <= 8; i += 4)
    {
        w[i + 8] = exp((-i) / temperature);
    }
}

void Ising::Metropolis()
{
    int deltaE;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            int rand_x = r.generateMersenneTwisterNumber(0.0, 1.0) * L;
            int rand_y = r.generateMersenneTwisterNumber(0.0, 1.0) * L;
            //std::cout << rand_x << ", ";
            //std::cout << rand_y << ", ";

            deltaE = getDeltaE(rand_x, rand_y);
            //std::cout << deltaE << " ";

            double random_check = r.generateMersenneTwisterNumber(0.0, 1.0);
            if (random_check <= w[deltaE + 8])
            {
                //std::cout << random_check << " : " << w[deltaE + 8] << " \n";
                spin_matrix[rand_x * L + rand_y] *= -1;

                M += 2 * spin_matrix[rand_x * L + rand_y];
                E += deltaE;
            }
        }
        //std::cout << "\n";
    }
    //std::cout << "\nEnergy = " << E << ", Magnetization = " << M << "\n";
}

int Ising::getDeltaE(int row, int col)
{
    int up = spin_matrix[row * L + periodic(col, 1)];
    int down = spin_matrix[row * L + periodic(col, -1)];
    int left = spin_matrix[periodic(row, -1) * L + col];
    int right = spin_matrix[periodic(row, 1) * L + col];
    return (2 * spin_matrix[row * L + col] * (up + down + left + right));
}

void Ising::MonteCarloCycler(int cycles)
{
    for (int cycle = 1; cycle < cycles + 1; cycle++)
    {
        Metropolis();

        expectation[0] += E;
        expectation[1] += E * E;
        expectation[2] += M;
        expectation[3] += M * M;
        expectation[4] += fabs(M);

        if (cycle % 1000 == 0)
        {
            computeValuesOfInterest(cycle);
            writeToFile(filename);
        }
    }
    //computeValuesOfInterest(cycles);
    //writeToFile(filename);
}

void Ising::computeValuesOfInterest(int cycle)
{
    double norm = 1 / (double)cycle;
    double expectation_norm[5];
    for (int i = 0; i < 5; i++)
    {
        expectation_norm[i] = expectation[i] * norm;
    }

    double E_var = (expectation_norm[1] - pow(expectation_norm[0], 2)) / (double)N;

    E_mean = expectation_norm[0] / (double)N;

    M_mean = expectation_norm[2] / (double)N;
    M_mean_abs = expectation_norm[4] / (double)N;

    specific_heat = E_var / pow(temperature, 2);

    susceptibility = (expectation_norm[3] - pow(expectation_norm[2], 2)) / ((double)N * temperature);
}

void Ising::writeToFile(std::string filename)
{
    std::ofstream ofile;

    ofile.open(filename + ".txt", std::ios::app);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    //ofile << std::setprecision(8) << current_cycle;
    ofile << std::setw(15) << std::setprecision(8) << temperature;
    ofile << std::setw(15) << std::setprecision(8) << E_mean;
    ofile << std::setw(15) << std::setprecision(8) << M_mean;
    ofile << std::setw(15) << std::setprecision(8) << specific_heat;
    ofile << std::setw(15) << std::setprecision(8) << susceptibility;
    ofile << std::setw(15) << std::setprecision(8) << M_mean_abs << "\n";
    ofile.close();
}