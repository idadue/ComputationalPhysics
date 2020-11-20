#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

class ising
{
public:
	void initialize(int, arma::mat&, double&, double&);
	double energy(arma::mat);
	arma::sword pbc(arma::sword, arma::sword);

	arma::vec transitionProb(double);
	arma::vec metropolis(int, int, double, std::string);
	void output(int, int, arma::vec, double, std::ofstream&);
};
