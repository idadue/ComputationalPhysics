#include "jacobisolver.h"

JacobiSolver::JacobiSolver(const unsigned int n, const float rho, const int potential, const float omega) : n{n}, h{rho / double(n)}
{
	a = -(1.0 / pow(h, 2));
	d = 2.0 / pow(h, 2);

	epsilon = 1.0e-8;
	if (potential == 1)
	{
		matrixInit(getPotential());
	}
	else if (potential == 2 && omega != 0)
	{
		matrixInit(getPotential(omega));
	}
	else
	{
		matrixInit(NULL);
	}
}

std::tuple<int, int, double> JacobiSolver::findMaxNonDiagElement(const arma::mat &V)
{
	double max = 0;
	std::tuple<int, int, double> max_ij;
	for (arma::uword i = 0; i < n; i++)
	{
		for (arma::uword j = i + 1; j < n; j++)
		{
			if (fabs(V(i, j)) > max)
			{
				max = fabs(V(i, j));
				max_ij = {i, j, max};
			}
		}
	}
	return max_ij;
}
/*Implementation of the rotation matrix as described by the Jacobi eigenvalue algorithm
Coulumns of R is filled with eigenvectors and B is rotated.*/
void JacobiSolver::rotate(unsigned int k, unsigned int l)
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

	for (arma::uword i = 0; i < n; i++)
	{
		if (i != k && i != l)
		{
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
	}
}

arma::mat JacobiSolver::sortR(const arma::mat &R, const arma::uvec &indSorted)
{
	arma::mat RSorted;
	RSorted.copy_size(R);
	for (arma::uword j = 0; j < R.n_cols; j++)
	{
		for (arma::uword i = 0; i < R.n_rows; i++)
		{
			RSorted(i, j) = R(i, indSorted(j));
		}
	}

	return RSorted;
}

/*Uses rotate() and findMaxNonDiagElements() to compute the Jacobi eigenvalue algorithm, with a multiplier such that 
we can run for enough iterations to reach the accuracy we need.*/
void JacobiSolver::solve(const int multiplier)
{
	int iterations = 0;
	int max_iterations = n * multiplier;
	auto maxoffdiag = findMaxNonDiagElement(B);

	while (pow(std::get<2>(maxoffdiag), 2) > epsilon && iterations < max_iterations)
	{
		rotate(std::get<0>(maxoffdiag), std::get<1>(maxoffdiag));
		maxoffdiag = findMaxNonDiagElement(B);
		iterations++;
	}
	arma::uvec indSorted = arma::sort_index(B.diag());

	R = sortR(R, indSorted);
	printf("End of loop at %d iterations \n", iterations);
}

double JacobiSolver::solveWithTime(const int multiplier)
{
	clock_t start, finish;
	start = clock();
	solve(multiplier);
	finish = clock();
	return (double(finish) - double(start)) / double(CLOCKS_PER_SEC);
}

std::tuple<arma::vec, arma::mat, double> JacobiSolver::armadilloSolution()
{
	clock_t start, finish;
	arma::vec eigval;
	arma::mat eigvec;
	start = clock();
	arma::eig_sym(eigval, eigvec, A, "std");
	finish = clock();
	double timing = (double(finish) - double(start)) / double(CLOCKS_PER_SEC);
	return std::tuple<arma::vec, arma::mat, double>(eigval, eigvec, timing);
}

arma::vec JacobiSolver::getPotential()
{
	arma::vec v(n);
	for (arma::uword i = 0; i < n; i++)
	{
		v(i) = pow((double(i) + 1.0) * h, 2);
	}
	return v;
}

arma::vec JacobiSolver::getPotential(const float omega)
{
	arma::vec v(n);
	for (arma::uword i = 0; i < n; i++)
	{
		double rho = (double(i) + 1.0) * h;
		v(i) = pow(omega, 2) * pow(rho, 2) + 1 / rho;
	}
	return v;
}

arma::mat JacobiSolver::getTestR()
{
	arma::mat R_copy = R;
	R_copy = R_copy * R_copy.t();
	for (arma::uword i = 0; i < n; i++)
	{
		R_copy(i, i) -= R_copy.diag()(i);
	}
	return R_copy;
}

double JacobiSolver::compareTrace()
{
	return abs(arma::sum(A.diag()) - arma::sum(B.diag()));
}

arma::mat JacobiSolver::getR()
{
	return R;
}

arma::mat JacobiSolver::getLambda()
{
	return arma::sort(B.diag());
}

arma::mat JacobiSolver::getAnalyticEigVec()
{
	arma::mat analyticalEigVec(n, n);
	for (arma::uword j = 0; j < n; j++)
	{
		for (arma::uword i = 1; i < n + 1; i++)
		{
			analyticalEigVec(i - 1, j) = sin(i * (double(j) + 1.0) * pi / (double(n + 1)));
		}
	}
	return analyticalEigVec;
}

arma::vec JacobiSolver::getAnalyticEigVal()
{
	arma::vec analyticEigVal(n);
	for (arma::uword j = 1; j < n + 1; j++)
	{
		analyticEigVal(j - 1) = d + 2 * a * std::cos(j * pi / (double(n + 1)));
	}
	return analyticEigVal;
}

void JacobiSolver::display(const arma::mat &v)
{
	std::cout.precision(4);
	std::cout.setf(std::ios::fixed);
	std::cout.width(9);
	v.raw_print(std::cout);
	printf("\n");
}