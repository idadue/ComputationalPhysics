#include "ising.h"
#include <iostream>

int main(int argc, char* argv[])
{
  int mcc = std::atoi(argv[1]);
  arma::sword n_spins = std::atoi(argv[2]);
  double temperature_start = std::atof(argv[3]);
  double temperature_end = std::atof(argv[4]);
  int n_temperature = std::atoi(argv[5]);
  int random_start = std::atoi(argv[6]);
  std::string filename = argv[7];
  int compute_mcc = std::atoi(argv[8]);

  ising isi;
  /*
  arma::mat spin_mat = arma::zeros<arma::mat>(n_spins, n_spins);
  double E = 0;
  double M = 0;
  isi.initialize(n_spins, spin_mat, E, M, random_start);
  spin_mat.print();
  */

  isi.mc_temp(n_spins, mcc, temperature_start, temperature_end, n_temperature, filename, random_start);

  if (compute_mcc == 1)
  {
    isi.mc_mcc(n_spins, mcc, temperature_start, temperature_end, n_temperature, filename, random_start);
  }
  //arma::vec exp = isi.metropolis(n_spins, mcc, temperature, filename, random_start);

  //std::cout << "L = " << n_spins << ", T = " << temperature << std::endl;
  return 0;
}
