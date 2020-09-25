#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "jacobisolver.h"
#include "filehandler.h"
#include <time.h>

//TODO Implement more tests?

/*
Contains test for:
	-Orthonormality of eigenvectors
	-Trace(A) = sum(lambda_i)
*/

TEST_CASE("Running the jacobi algorithm") {
	unsigned int n, multiplier, potensial;
	float rho, omega;
	std::ifstream inData;
	inData.open("inData.txt");
	inData.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	inData >> n >> rho >> multiplier >> potensial >> omega;

	JacobiSolver j(n, rho, multiplier, potensial, omega);

	auto armaeig = j.armadilloSolution();
	auto R = j.getTestR();

	REQUIRE(abs(R.max()) < 1.0e-14);
	double trace_diff = j.compareTrace();
	REQUIRE(trace_diff < 1.0e-02);

	//printf("Analytical eigenvectors: \n");
	//std::cout << j.getAnalyticEigVec() << std::endl;
	//printf("Numerical eigenvectors(unsorted): \n");
	//std::cout << j.getR() << std::endl;
	//printf("Numerical eigenvectors(armadillo 'std'): \n");
	//std::cout << armaeig.second << std::endl;
	//printf("Test of orthonormality: \n");
	//std::cout << j.getR().t() * j.getR() << std::endl;

	//printf("Analytical eigenvalues: \n");
	//std::cout << j.getAnalyticEigVal() << std::endl;
	printf("Numerical eigenvalues(unsorted): \n");
	//std::cout << arma::sort(j.getLambda()) << std::endl;
	std::cout << j.getLambda() << std::endl;
	printf("Numerical eigenvalues(armadillo 'std'): \n");
	std::cout << armaeig.first << std::endl;

	//FileWriter file("eigenvalues");
	//file.writeToFile(j.getLambda());
}