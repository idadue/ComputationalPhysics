#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "jacobisolver.h"
#include "filehandler.h"

/*
Contains test for:
	-Orthonormality of eigenvectors
	-Trace(A) = sum(B.diag())
*/

TEST_CASE("Running the jacobi algorithm")
{
	unsigned int n, multiplier, potential;
	float rho, omega;
	std::string folder, doArma;
	std::ifstream inData;
	int counter = 0;
	std::vector<double> timing;

	inData.open("inData.txt");
	inData.ignore(std::numeric_limits<std::streamsize>::max(), '\\');

	while (!inData.eof())
	{
		inData >> n >> rho >> multiplier >> potential >> omega >> folder >> doArma;
		std::cout << "Using n = " << n << ", rho = " << rho << ", potential = " << potential << ", and  omega = " << omega << "\n";

		FileHandler file(folder);
		std::string filename = std::to_string(n);
		std::string filenameEigVec = std::to_string(n) + "_eigenvector1";
		std::string filenameArma = "arma_" + std::to_string(n);

		JacobiSolver j(n, rho, potential, omega);
		double time = j.solveWithTime(multiplier);

		arma::mat testR = j.getTestR();
		REQUIRE(abs(testR.max()) < 1.0e-14);
		double trace_diff = j.compareTrace();
		REQUIRE(trace_diff < 1.0e-02);

		if (folder != ".")
		{
			file.writeRawToFile(filename, j.getLambda());
			file.writeToFile(filenameEigVec, j.getR().col(0));

			if (doArma == "arma")
			{
				std::tuple<arma::vec, arma::mat, double> armaeig = j.armadilloSolution();
				file.writeRawToFile(filenameArma, std::get<0>(armaeig));
				std::cout << "Execution time of arma algorithm is: " << std::get<2>(armaeig) << " seconds \n";
			}

			if (potential == 0)
			{
				arma::vec analytic_eigval = j.getAnalyticEigVal();
				arma::vec analytic_eigvec = j.getAnalyticEigVec().col(0);

				std::string a_eigval = "analytic_eigval_" + std::to_string(n);
				std::string a_eigvec = "analytic_eigvec1_" + std::to_string(n);

				file.writeToFile(a_eigval, analytic_eigval);
				file.writeToFile(a_eigvec, analytic_eigvec);
			}
		}
		else
		{
			timing.push_back(time);
			counter++;
		}

		std::cout << "Execution time of algorithm is: " << time << " seconds \n \n";
	}

	if (!timing.empty())
	{
		arma::vec t = timing;
		FileHandler file("bbeam");
		std::string file_t = "timing";
		file.writeToFile(file_t, timing);
		printf("Average execution time for buckling beam problem with n = 100 was t = %f seconds \n", ((arma::sum(t)) / double(counter)));
	}

	inData.close();
}