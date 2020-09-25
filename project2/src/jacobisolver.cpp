#include "jacobisolver.h"

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

*/

JacobiSolver::JacobiSolver(const unsigned int n, const float rho, const int potential, const float omega) : n{ n }, h{ rho / double(n) }
{
	a = -(1.0 / pow(h, 2));
	d = 2.0 / pow(h, 2);

	epsilon = 1.0e-12;
	if (potential == 1) {
		matrixInit(getPotential());
	}
	else if (potential == 2 && omega != 0) {
		matrixInit(getPotential(omega));
	}
	else {
		matrixInit(NULL);
	}
}

JacobiSolver::JacobiSolver(const unsigned int n, const float rho, const int multiplier, const int potential, const float omega) : JacobiSolver(n, rho, potential, omega)
{
	solve(multiplier);
}

std::tuple<int, int, double> JacobiSolver::findMaxNonDiagElement(const arma::mat& V)
{
	double max = 0;
	std::tuple<int, int, double> max_ij;
	for (arma::uword i = 0; i < n; i++) {
		for (arma::uword j = 0; j < n; j++) {
			if (i != j && fabs(V(i, j)) > max) {
				max = fabs(V(i, j));
				max_ij = { i,j, max };
			}
		}
	}
	return max_ij;
}

void JacobiSolver::rotate(int k, int l)
{
	double s, c;

	if (B(k, l) != 0.0)
	{
		double t, tau;
		tau = (B(l, l) - B(k, k)) / (2 * B(k, l));

		(tau >= 0.0) ? t = 1.0 / (tau + sqrt(1.0 + pow(tau, 2))) : t = 1.0 / (tau - sqrt(1.0 + pow(tau, 2)));

		c = 1.0 / sqrt(1 + pow(t, 2));
		s = c * t;
	}
	else
	{
		c = 1.0;
		s = 0.0;
	}
	double a_ik, a_il, r_ik, r_il;
	double a_kk = B(k, k);
	double a_ll = B(l, l);
	B(k, k) = c * c * a_kk - 2.0 * c * s * B(k, l) + s * s * a_ll;
	B(l, l) = s * s * a_kk + 2.0 * c * s * B(k, l) + c * c * a_ll;
	B(k, l) = 0;
	B(l, k) = 0;

	for (unsigned int i = 0; i < n; i++) {
		if (i != k && i != l) {
			a_ik = B(i, k);
			a_il = B(i, l);
			B(i, k) = c * a_ik - s * a_il;
			B(k, i) = B(i, k);
			B(i, l) = c * a_il + s * a_ik;
			B(l, i) = B(i, l);
		}
		r_ik = R(i, k);
		r_il = R(i, l);
		R(i, k) = c * r_ik - s * r_il;
		R(i, l) = c * r_il + s * r_ik;
		lambda(i) = B(i, i);
	}
}

void JacobiSolver::solve(int multiplier)
{
	int iterations = 0;
	int max_iterations = n * multiplier;
	auto maxoffdiag = findMaxNonDiagElement(B);

	while (pow(std::get<2>(maxoffdiag), 2) > epsilon && iterations < max_iterations) {
		rotate(std::get<0>(maxoffdiag), std::get<1>(maxoffdiag));
		maxoffdiag = findMaxNonDiagElement(B);
		iterations++;
	}
	std::cout << "End of loop at " << iterations << " iterations completed.\n \n";
	arma::normalise(R, 1);
}

arma::vec JacobiSolver::getPotential()
{
	arma::vec v(n);
	for (arma::uword i = 0; i < n; i++) {
		v(i) = pow((double(i) + 1.0) * h, 2);
	}
	return v;
}

arma::vec JacobiSolver::getPotential(const float omega)
{
	arma::vec v(n);
	for (arma::uword i = 0; i < n; i++) {
		double rho = (double(i) + 1.0) * h;
		v(i) = pow(omega, 2) * pow(rho, 2) + 1 / rho;
	}
	return v;
}

std::pair<arma::vec, arma::mat> JacobiSolver::armadilloSolution()
{
	arma::vec eigval;
	arma::mat eigvec;
	arma::eig_sym(eigval, eigvec, A, "std");
	return std::pair<arma::vec, arma::mat>(eigval, eigvec);
}

arma::mat JacobiSolver::getTestR()
{
	arma::mat R_copy = R;
	R_copy = R_copy * R_copy.t();
	for (unsigned int i = 0; i < n; i++) {
		R_copy(i, i) -= R_copy.diag()(i);
	}
	return R_copy;
}

double JacobiSolver::compareTrace()
{
	double A_diag = arma::sum(A.diag());
	return abs(A_diag - arma::sum(lambda));
}

arma::mat JacobiSolver::getA()
{
	return A;
}

arma::mat JacobiSolver::getR()
{
	return R;
}

arma::mat JacobiSolver::getLambda()
{
	return lambda;
}

arma::mat JacobiSolver::getAnalyticEigVec()
{
	arma::mat analyticalEigVec(n, n);
	for (arma::uword j = 0; j < n; j++) {
		for (arma::uword i = 1; i < n + 1; i++) {
			analyticalEigVec(i - 1, j) = sin(i * (double(j) + 1.0) * pi / (double(n)));
		}
	}
	//analyticalEigVec = arma::normalise(analyticalEigVec, 1);
	return analyticalEigVec;
}

arma::mat JacobiSolver::getAnalyticEigVal()
{
	arma::vec analyticEigVal(n);
	for (unsigned int j = 1; j < n + 1; j++)
	{
		analyticEigVal(j - 1) = d + 2 * a * std::cos(j * pi / (double(n)));
	}
	return analyticEigVal;
}