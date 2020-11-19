#include <armadillo>
#include <cmath>
#include <cstdlib>

class ising
{
public:
	void initialize(int, arma::mat &, double&, double&);
	double energy(arma::mat);
	arma::sword pbc(arma::sword, arma::sword);

	arma::vec transitionProb(double);
	arma::vec metropolis(int, int, double);
	void output(int, arma::vec, double);
};
