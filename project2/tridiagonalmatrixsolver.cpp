#include "TridiagonalMatrixSolver.h"

double TridiagonalMatrixSolver::f(double x)
{
	return 100 * exp(-10 * x);
}

void TridiagonalMatrixSolver::display()
{
	std::cout.precision(7);
	std::cout.setf(std::ios::fixed);
	v.raw_print(std::cout, "v:");
}

void TridiagonalMatrixSolver::writeToFile(std::string filename)
{
	outputFile.open(filename, std::fstream::out);
	if (outputFile.is_open()) {
		outputFile << 0 << "\n";
		for (arma::uword i = 0; i < n; i++) {
			outputFile << v[i] << "\n";
		}
		outputFile << 0;
	}
	outputFile.close();
}

/*------------------------------------------------------------*/

void ThomasSolver::forward_substitution()
{
	for (arma::uword i = 1; i < n; i++) {
		d[i] += -c * a / d[i - 1];
		b_tilde[i] += -b_tilde[i - 1] * a / d[i - 1];
	}
}

void ThomasSolver::backward_substitution()
{
	v[n - 1] = b_tilde[n - 1] / d[n - 1];
	for (int i = n - 2; i > -1; i--) {
		v[i] = (b_tilde[i] - c * v[i + 1]) / d[i];
	}
}

double ThomasSolver::getExecutionTime()
{
	return (finish - start) / CLOCKS_PER_SEC;
}

void ThomasSolver::solve()
{
	start = clock();
	forward_substitution();
	backward_substitution();
	finish = clock();
}

/*------------------------------------------------------------*/

void SpecialThomasSolver::forward_subtitution()
{
	for (arma::uword i = 0; i < n; i++) {
		//less optimal? probably...
		d[i] = (i + 2.0) / (i + 1.0);
		b_tilde[i] = (i >= 1) ? (step_size_sqrd * f((i + 1.0) * step_size) + b_tilde[i - 1] / d[i - 1]) :
			step_size_sqrd * f(step_size);
	}
}

void SpecialThomasSolver::backward_subtitution()
{
	v[n - 1] = b_tilde[n - 1] / d[n - 1];
	for (int i = n - 2; i > -1; i--) {
		v[i] = (b_tilde[i] + v[i + 1]) / d[i];
	}
}

void SpecialThomasSolver::solve()
{
	forward_subtitution();
	backward_subtitution();
}

/*------------------------------------------------------------*/

void LUSolver::solve()
{
	for (arma::uword i = 0; i < n; i++) {
		v[i] = step_size_sqrd * f((i + 1.0) * step_size);
		for (arma::uword j = 0; j < i; j++) {
			v[i] -= L(i, j) * v[j];
		}
	}

	for (int i = n - 1; i > -1; i--) {
		for (int j = n - 1; j > i; j--) {
			v[i] -= U(i, j) * v[j];
		}
		v[i] /= U(i, i);
	}
}