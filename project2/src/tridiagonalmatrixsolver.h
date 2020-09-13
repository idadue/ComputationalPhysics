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

class TridiagonalMatrixSolver
{
protected:
	unsigned int n;
	double step_size, step_size_sqrd;
	clock_t start, finish;
	arma::mat v;
	std::fstream outputFile;
public:
	TridiagonalMatrixSolver(const unsigned int n)
	{
		this->n = n;

		v.zeros(n);

		step_size = 1.0 / (double(n) + 1);
		step_size_sqrd = step_size * step_size;
	}

	double f(double x);
	virtual void display();
	virtual void writeToFile(std::string filename);
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
	ThomasSolver(const unsigned int n, const double a, const double b, const double c) : TridiagonalMatrixSolver(n) {
		this->d = this->d.ones(n) * b;
		this->a = a; this->c = c;

		b_tilde.zeros(n);

		for (arma::uword i = 0; i < n; i++) {
			b_tilde[i] = step_size_sqrd * f((i + 1.0) * step_size);
		}
	}
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
	SpecialThomasSolver(const unsigned int n) : TridiagonalMatrixSolver(n) {
		d = d.ones(n) * 2;
		b_tilde.zeros(n);
	}

	void solve();
};

/*------------------------------------------------------------*/

class LUSolver : public TridiagonalMatrixSolver {
private:
	arma::mat A, L, U;
public:
	LUSolver(const unsigned int n, const int a, const int b, const int c) : TridiagonalMatrixSolver(n) {
		A.zeros(n, n);
		A.diag() += b;
		A.diag(-1) += a;
		A.diag(1) += c;

		arma::lu(L, U, A);
	}

	void solve();
};

#endif