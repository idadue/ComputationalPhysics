#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

class ising
{
public:
	void initialize(int, arma::mat&, double&, double&, int);
	arma::ivec pbc;
	double compute_energy(arma::mat);

	arma::vec transitionProb(double);

	void mc_L(int, double, double, int, std::string, double, int, int, int, int);
	void mc_temp(int, int, double, double, int, std::string, double, int);
	void mc_mcc(int, int, double, double, int, std::string, double, int);
	void mc_config_temp(int, int, double, double, int, std::string, double, int);
	std::ofstream out_thermo;
	std::ofstream out_energy;
	std::ofstream out_temp;
	std::ofstream out_mcc;
	int flip_count;

	arma::vec metropolis(int, int, double, std::string, double, int);
	arma::vec data_conversion(int, int, arma::vec, double);
	void output(arma::vec, std::ofstream&);
};
