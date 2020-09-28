#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include "jacobisolver.h"
#include "filehandler.h"
#include <time.h>

//TODO Implement more tests?

/*
Contains test for:
	-Orthonormality of eigenvectors
	-Trace(A) = sum(B.diag())
*/

TEST_CASE("Running the jacobi algorithm") {
	unsigned int n, multiplier, potential;
	float rho, omega;
	std::string folder;
	std::ifstream inData;
	int counter = 0;
	std::vector<double> timing;

	inData.open("inData.txt");
	inData.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	inData.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	while (!inData.eof()) {
		inData >> n >> rho >> multiplier >> potential >> omega >> folder;

		FileHandler file(folder);

		std::string filename = std::to_string(n);
		std::string filenameEigVec = std::to_string(n) + "_eigenvector1";
		std::string filenameArma = "arma_" + std::to_string(n);

		JacobiSolver j(n, rho, potential, omega);
		double time = j.solveWithTime(multiplier);

		std::pair<arma::vec, arma::mat> armaeig = j.armadilloSolution();
		arma::mat testR = j.getTestR();
		REQUIRE(abs(testR.max()) < 1.0e-14);
		double trace_diff = j.compareTrace();
		REQUIRE(trace_diff < 1.0e-02);

		if (folder != ".") {
			file.writeToFile(filename, j.getLambda());
			file.writeToFile(filenameEigVec, j.getR().col(0));
			file.writeToFile(filenameArma, armaeig.first);

			if (potential == 0) {
				arma::vec analytic_eigval = j.getAnalyticEigVal();
				arma::vec analytic_eigvec = j.getAnalyticEigVec().col(0);

				std::string a_eigval = "analytic_eigval_" + std::to_string(n);
				std::string a_eigvec = "analytic_eigvec1_" + std::to_string(n);

				file.writeToFile(a_eigval, analytic_eigval);
				file.writeToFile(a_eigvec, analytic_eigvec);
			}
		}
		else {
			timing.push_back(time);
			counter++;
		}

		std::cout << "Using n = " << n << ", rho = " << rho << ", potential = " <<
			potential << ", and  omega = " << omega << "\n";
		std::cout << "Execution time of algorithm is: " << time << " seconds \n \n";
	}
	arma::vec t = timing;
	if (counter != 0) {
		printf("Average execution time for buckling beam problem with n = 100 was t = %f seconds \n", ((arma::sum(t)) / double(counter)));
	}

	FileHandler file("bbeam");
	std::string file_t = "timing";
	file.writeToFile(file_t, timing);

	/*
	printf("Analytical eigenvectors: \n");
	j.display(j.getAnalyticEigVec());
	printf("Numerical eigenvectors(unsorted): \n");
	j.display(j.getR());
	printf("Numerical eigenvectors(armadillo 'std'): \n");
	j.display(armaeig.second);*/
	//printf("Test of orthonormality: \n");
	//std::cout << j.getR().t() * j.getR() << std::endl;

	//printf("Analytical eigenvalues: \n");
	//j.display(j.getAnalyticEigVal());
	//printf("Numerical eigenvalues(unsorted): \n");
	//std::cout << arma::sort(j.getLambda()) << std::endl;
	//j.display(j.getLambda());
	//printf("Numerical eigenvalues(armadillo 'std'): \n");
	//j.display(armaeig.first);

	inData.close();
}