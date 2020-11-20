#include "ising.h"
#include <iostream>

int main(int argc, char* argv[])
{
  int mcc = std::atoi(argv[1]);
  arma::sword n_spins = std::atoi(argv[2]);
  double temperature = std::atof(argv[3]);
  std::string filename = argv[4];

  ising isi;

  //arma::mat spin_mat = arma::zeros<arma::mat>(n_spins, n_spins);
  //double E = 0;
  //double M = 0;

  //std::cout << spin_mat << std::endl;
  //isi.initialize(n_spins, spin_mat, E, M);
  //std::cout << spin_mat << std::endl;


  //isi.transitionProb(temperature).print();

  arma::vec exp = isi.metropolis(n_spins, mcc, temperature, filename);
  return 0;
}
