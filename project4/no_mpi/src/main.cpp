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
  int n_L = std::atoi(argv[9]);
  int spin_end = std::atoi(argv[10]);
  int equilibrium_factor = std::atof(argv[11]);

  ising isi;

  if (compute_mcc == 1)
  {
    isi.mc_mcc(n_spins, mcc, temperature_start, temperature_end, n_temperature, filename, equilibrium_factor, random_start);
  }
  else if (compute_mcc == 2)
  {
    isi.mc_config_temp(n_spins, mcc, temperature_start, temperature_end, n_temperature, filename, equilibrium_factor, random_start);
  }
  else if (n_L > 1)
  {
    isi.mc_L(mcc, temperature_start, temperature_end, n_temperature, filename, equilibrium_factor, random_start, n_spins, n_L, spin_end);
  }
  else
  {
    isi.mc_temp(n_spins, mcc, temperature_start, temperature_end, n_temperature, filename, equilibrium_factor, random_start);
  }
  return 0;
}
