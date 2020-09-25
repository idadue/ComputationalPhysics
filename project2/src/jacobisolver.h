#pragma once

#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>

//TODO: Remove surplous methods, variables
const long double pi = 3.141592653589793238462643383279502884;

class JacobiSolver
{
private:
	unsigned int n;
	double a, d;
	double epsilon;
	double h;
	float rho;
	arma::mat A, B, R;
	arma::vec lambda;

	std::tuple<int, int, double> findMaxNonDiagElement(const arma::mat& V);

	void rotate(int k, int l);

public:
	JacobiSolver(const unsigned int n, const float rho, const int potential, const float omega);
	JacobiSolver(const unsigned int n, const float rho, const int multiplier, const int potential, const float omega);

	template <class T>
	void matrixInit(T potential);

	void solve(int multiplier);

	arma::vec getPotential();
	arma::vec getPotential(const float omega);

	std::pair<arma::vec, arma::mat> armadilloSolution();

	//Test
	arma::mat getTestR();
	double compareTrace();

	arma::mat getA();
	arma::mat getR();
	arma::mat getAnalyticEigVec();
	arma::mat getAnalyticEigVal();
	arma::mat getLambda();
};

template<class T>
inline void JacobiSolver::matrixInit(T potential)
{
	A = A.zeros(n, n);
	(typeid(potential) == typeid(arma::vec)) ? A.diag() += potential + d : A.diag() += d;
	A.diag(-1) += a;
	A.diag(1) += a;

	B = A;

	R = R.zeros(n, n);
	R.diag() += 1;

	lambda = lambda.zeros(n);
}
