#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

class disease
{
public:
  void initialize(double, int, int, int, double, double, double, double, double, double, double, double, double, double);
  arma::vec SIR = arma::zeros<arma::vec>(6);
  arma::vec EXP = arma::zeros<arma::vec>(6);

  int mcc;

  int total;
  double a_0;
  double a;
  double b;
  double c;
  double d;
  double d_I;
  double e;
  double dt;
  double time;
  double f;

  void transitionProb();
  // Transition probability from SIRDDU to SIRDDU
  arma::mat tProb = arma::zeros<arma::mat>(6, 6);

  void transmissionRate(double);
  double A;
  double w;
  void vaccinationRate(double);
  double f_0;
  double vaccineTime;

  void monteCarlo();
  std::ofstream out_exp;
  std::ofstream out_val;
  void output(arma::vec, std::ofstream&);
  void solver(std::string);
};
