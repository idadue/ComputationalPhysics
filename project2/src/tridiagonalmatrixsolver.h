#ifndef TriDiagonalMatrixSolver_h
#define TriDiagonalMatrixSolver_h

//Inspired by https://github.com/reneaas/classes_in_cpp/tree/master/advanced_tutorial

#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include <time.h>

#define pi std::_Pi

//Note: [] on armadillo objects != (). Safe for use on one dimensional objects, but NOT multi dimensional objects
//usigned int on backwards substitution loops is not good
//Member initialization list does not work with armadillo vec class. setting #define ARMA_USE_CXX11 in config.hpp should fix this, although testing
//shows this is at best unstable. Recommend avoiding or using mat instead, although mat does not seem to work when used as a matrix

class TridiagonalMatrixSolver
{
protected:
	unsigned int n;
	double step_size, step_size_sqrd;
	clock_t start, finish;
	arma::mat v; // d, b_tilde;
	std::fstream outputFile;
public:
	TridiagonalMatrixSolver(const unsigned int n);
	//TridiagonalMatrixSolver(const unsigned int n, const unsigned a);

	double f(double x);
	virtual void display();
	virtual void writeToFile(std::string filename);

	//Add your code here, repeat for other classes.
};

/*------------------------------------------------------------*/

class ThomasSolver : public TridiagonalMatrixSolver
{
private:
	arma::vec d, b_tilde;
	double a, c;
	void forward_substitution();
	void backward_substitution();
public:
	ThomasSolver(const unsigned int n, const double a, const double b, const double c);
	double getExecutionTime();
	void solve();
};

/*------------------------------------------------------------*/

class SpecialThomasSolver : public TridiagonalMatrixSolver {
private:
	arma::vec d, b_tilde;

	void forward_subtitution();
	void backward_subtitution();
public:
	SpecialThomasSolver(const unsigned int n);

	void solve();
};

/*------------------------------------------------------------*/

class LUSolver : public TridiagonalMatrixSolver {
private:
	arma::mat A, L, U;
public:
	LUSolver(const unsigned int n, const int a, const int b, const int c);

	void solve();
};

/*------------------------------------------------------------*/

class JacobiSolver : TridiagonalMatrixSolver {
private:
	//will probably change some of these names and function types
	double* analyticalEigenval;
	//vector of vectors?
	double* analyticalEigenvec;
	double a, d;
	double step_size, step_size_sqrd;
	unsigned int rho;
	arma::mat A;
public:
	JacobiSolver(const int n, const float rho);
	~JacobiSolver();
	void analyticalEigenvalues();
	//not tested
	void analyticalEigenvector();

	double* getanalyticalEigenvalues();
};

#endif