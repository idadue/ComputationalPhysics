#include "ising.h"

//initialize energy, magnetization
void ising::initialize(int n_spins, arma::mat & spin_mat, double& E, double& M)
{
  //setup spin matrix
  for (arma::sword i = 0; i < n_spins; i++){
    for (arma::sword j = 0; j < n_spins; j++){
      spin_mat(i,j) = 1;
      M += (float) spin_mat(i,j);
    }
  }
  E = energy(spin_mat);
}

// Calculates the energy of a lattice
double ising::energy(arma::mat spin_mat)
{
  double energy = 0;
  arma::sword n_spins = spin_mat.n_rows;
  for (arma::sword i = 0; i < n_spins; i++){
      for (arma::sword j = 0; j < n_spins; j++){
        energy -= spin_mat(i,j) * spin_mat(pbc(i + 1, n_spins), j) + spin_mat(i, j) * spin_mat(i, pbc(j + 1, n_spins));
      }
    }
  return energy;
}

// implement periodic boundary conditions
arma::sword ising::pbc(arma::sword coordinate, arma::sword max_coordinate)
{
  if (coordinate >= max_coordinate)
  {
    coordinate = coordinate - max_coordinate;
  }
  else if (coordinate < 0)
  {
    coordinate = coordinate + max_coordinate;
  }
  return coordinate;
}

// Returns the transition probability for positive energy differences
arma::vec ising::transitionProb(double temperature)
{
  int n_ediff = 5;
  arma::vec probability = arma::zeros<arma::vec>(n_ediff);
  for (arma::sword i = 0; i < 5; i++)
  {
    //std::cout << -4*(i-2.0) << std::endl;
    probability(i) = std::exp(-4*(i-2.0)/temperature);
  }
  return probability;
}

arma::vec ising::metropolis(int n_spins, int mcc, double temperature)
{
  //initializing values
  double energy = 0.;
  double dEnergy;
  double mag_mom = 0.;
  arma::vec exp = arma::zeros<arma::vec>(5);
  arma::mat spin_mat = arma::zeros<arma::mat>(n_spins, n_spins);
  initialize(n_spins, spin_mat, energy, mag_mom);

  // Initialize rng
  std::random_device rd;
  std::mt19937_64 gen(rd());

  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

  // Set up transition probability vector
  arma::vec probVec = transitionProb(temperature);
  probVec.print();
  for (int cycles = 1; cycles < mcc; cycles++)
  {
    std::cout << "Cycle #" << cycles << std::endl;
    for (arma::uword i = 0; i < n_spins; i++)
    {
      for (arma::uword j = 0; j < n_spins; j++)
      {
        //std::cout << "Coordinates " << i << ", " << j << std::endl;
        arma::sword i_flip = (arma::sword) (RandomNumberGenerator(gen)*(double)n_spins);
        arma::sword j_flip = (arma::sword) (RandomNumberGenerator(gen)*(double)n_spins);

        //std::cout << i_flip << ", " << j_flip << std::endl;

        dEnergy = 2 * spin_mat(i_flip, j_flip) * (spin_mat(pbc(i_flip + 1, n_spins), j_flip) + spin_mat(pbc(i_flip - 1, n_spins), j_flip) + spin_mat(i_flip, pbc(j_flip + 1, n_spins)) + spin_mat(i_flip, pbc(j_flip - 1, n_spins)));

        arma::sword energyCoordinate = (arma::sword) dEnergy/4 + 2;
        //std::cout << "denergy = " << dEnergy << std::endl;
        //std::cout << "energy coordinate = " << energyCoordinate << std::endl;
        //std::cout << "transitionProb = " << probVec(energyCoordinate) << std::endl;

        if (RandomNumberGenerator(gen) < probVec(energyCoordinate))
        {
          spin_mat(i_flip, j_flip) *= -1;
          energy += dEnergy;
          mag_mom += 2*spin_mat(i_flip, j_flip);
          //std::cout << "Flipped " << i_flip << ", " << j_flip << std::endl;
        }
      }
    }
    if (cycles > mcc*0.1)
    {
      exp(0) += energy;
      exp(1) += energy*energy;
      exp(2) += mag_mom;
      exp(3) += mag_mom*mag_mom;
      exp(4) += std::abs(mag_mom);
    }
  }
  exp /= (mcc*0.9);
  std::cout << "Fission Mailed Successfully" << std::endl;
  return exp;
}

void ising::output(int n_spins, arma::vec exp, double temperature)
{
  double exp_energy = exp(0);
  double var_energy = (exp(1) - exp(0)*exp(0));
  double exp_magmom = exp(2);
  double var_magmom = (exp(3) - exp(2)*exp(2));

  double spec_heat = var_energy / temperature / temperature;
  double mag_suscept = var_magmom / temperature;
  std::cout << "E(E) = " << exp_energy << ", Var(E) = " << var_energy << std::endl;
  std::cout << "E(M) = " << exp_magmom << ", Var(M) = " << var_magmom << std::endl;
  std::cout << "C_V = " << spec_heat << ", X = " << mag_suscept << std::endl;

}
