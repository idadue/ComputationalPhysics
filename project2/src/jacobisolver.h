#pragma once
#include <armadillo>
#include <cmath>
#include <iostream>
#include <iomanip>

/*
Class for solving eigenvalue problems using the Jacobi eigenvalue algorithm.

Variables:
	h - step size
	h_sqrd - h^2

	A - Tridiagonal matrix (Toeplitz)
	B - Copy of A, used in rotation
	R - Matrix containing calculated eigenvectors

	a, d - problem specific constants
	epsilon - tolerance
	potential - int, choose what hard coded potential to use
	omega - frequency, only used with the corresponding potential

Functions:
	void matrixInit() - takes any type T, altough in reality only type arma::vec is accounted for.
	std::tuple findMaxNonDiagElement - returns the biggest absoulute non diagonal element of a matrix with indices.

	void solve() - Solve for eigenvalues of tridiagonal Toeplitz matrix using Jacobi eigenvalue algorithm
	double solveWithTime() - Same as solve() but returns the time it took to completed computation

	void rotate() - the actual method for implementing the Jacobi eigenvalue algorithm
	arma::mat sortR - sorts the eigenvector matrix to match with the sorted eigenvalues

	arma::vec getPotential() - returns a potential function to be added to the super and sub diagonals elements of the Toeplitz matrix.
	arma::mat getTestR() - Returns an R matrix with diagonals removed so the orthonormality of R can be tested.
	double compareTrace() - Compares the trace of A with the sum of the diagonal elements of B;

	Other functions should be self explanatory.

*/

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

	std::tuple<int, int, double> findMaxNonDiagElement(const arma::mat& V);

	void rotate(int k, int l);
	arma::mat sortR(const arma::mat& R, const arma::uvec& indSorted);

public:
	JacobiSolver(const unsigned int n, const float rho, const int potential, const float omega);

	template <class T>
	void matrixInit(T potential);

	void solve(const int multiplier);
	double solveWithTime(const int multiplier);

	std::pair<arma::vec, arma::mat> armadilloSolution();

	arma::vec getPotential();
	arma::vec getPotential(const float omega);

	//Test
	arma::mat getTestR();
	double compareTrace();

	arma::mat getR();
	arma::mat getLambda();
	arma::mat getAnalyticEigVec();
	arma::vec getAnalyticEigVal();

	void display(const arma::mat& v);
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
}
