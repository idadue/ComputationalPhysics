#include "ising.h"
#include <iostream>

int main(int argc, char* argv[])
{
  ising isi;
  arma::sword n_spins = 2;
  arma::mat spin_mat = arma::zeros<arma::mat>(n_spins, n_spins);
  double E = 0;
  double M = 0;

  //std::cout << spin_mat << std::endl;
  isi.initialize(n_spins, spin_mat, E, M);
  //std::cout << spin_mat << std::endl;

  int mcc = std::atoi(argv[1]);
  double temperature = 1;

  isi.transitionProb(temperature).print();
  arma::vec exp = isi.metropolis(n_spins, mcc, temperature);
  exp.print();
  isi.output(n_spins, exp, temperature);
  return 0;
}
