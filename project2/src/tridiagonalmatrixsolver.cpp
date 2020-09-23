#include "TridiagonalMatrixSolver.h"

TridiagonalMatrixSolver::TridiagonalMatrixSolver(const unsigned int n) : n{ n }, start{ 0 }, finish{ 0 }, v{ v.zeros(n) }
{
	//initializing arma::mat should probably be done directly in constructor
	//v = v.zeros(n);
	step_size = 1.0 / (double(n) + 1);
	step_size_sqrd = step_size * step_size;
}

double TridiagonalMatrixSolver::f(double x)
{
	return 100 * exp(-10 * x);
}

void TridiagonalMatrixSolver::display()
{
	std::cout.precision(7);
	std::cout.setf(std::ios::fixed);
	v.raw_print(std::cout, "v:");
	printf("\n");
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

ThomasSolver::ThomasSolver(const unsigned int n, const double a, const double b, const double c) : TridiagonalMatrixSolver(n), a{ a }, c{ c }
{
	b_tilde.zeros(n);
	d = d.ones(n) * b;
	for (arma::uword i = 0; i < n; i++) {
		b_tilde[i] = step_size_sqrd * f((i + 1.0) * step_size);
	}
}

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

SpecialThomasSolver::SpecialThomasSolver(const unsigned int n) : TridiagonalMatrixSolver(n) {
	d = d.ones(n) * 2;
	b_tilde.zeros(n);
}

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

LUSolver::LUSolver(const unsigned int n, const int a, const int b, const int c) : TridiagonalMatrixSolver(n)
{
	A = A.zeros(n, n);
	A.diag() += b;
	A.diag(-1) += a;
	A.diag(1) += c;

	arma::lu(L, U, A);
}

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

/*------------------------------------------------------------*/

JacobiSolver::JacobiSolver(const int n, const float rho) : TridiagonalMatrixSolver(n)
{
	step_size = rho / double(n);
	step_size_sqrd = step_size * step_size;

	a = -(1 / step_size_sqrd);
	d = 2 / step_size_sqrd;

	A = A.zeros(n, n);
	A.diag() += d;
	A.diag(-1) += a;
	A.diag(1) += a;
}

JacobiSolver::~JacobiSolver() {
	delete[] analyticalEigenval;
	delete[] analyticalEigenvec;
}

void JacobiSolver::analyticalEigenvalues()
{
	analyticalEigenval = new double[n];
	for (unsigned int j = 1; j < n + 1; j++) {
		analyticalEigenval[j - 1] = d + 2 * a * std::cos(j * pi / (double(n) + 1));
	}
}

void JacobiSolver::analyticalEigenvector()
{
	analyticalEigenvec = new double(n);
	for (unsigned int j = 1; j < n + 1; j++) {
		analyticalEigenvec[j - 1] = std::sin(j * pi / (double(n) + 1));
	}
}

double* JacobiSolver::getanalyticalEigenvalues()
{
	return analyticalEigenval;
}