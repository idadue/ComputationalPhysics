#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

class ising
{
public:
	void initialize(int, arma::mat&, double&, double&, int);
	double compute_energy(arma::mat);
	arma::sword pbc(arma::sword, arma::sword);

	arma::vec transitionProb(double);

	void mc_temp(int, int, double, double, int, std::string, int);
	void mc_mcc(int, int, double, double, int, std::string, int);
	std::ofstream out_thermo;
	std::ofstream out_energy;
	std::ofstream out_temp;
	std::ofstream out_mcc;
	int flip_count;

	arma::vec metropolis(int, int, double, std::string, int);
	bool is_equilibrium();
	arma::vec data_conversion(int, int, arma::vec, double);
	void output(arma::vec, std::ofstream&);
};
