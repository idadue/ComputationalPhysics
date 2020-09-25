#pragma once
#include <fstream>
#include <string>
#include <armadillo>

class FileHandler {
private:
	std::ofstream outFile;
	std::ifstream inFile;
	std::string filename;

public:
	FileHandler() : filename{ "default.csv" } {}
	FileHandler(std::string filename) : filename{ filename + ".csv" } {}

	int writeMatrixToFile(const arma::mat& V) {
		fileOpener(exists(filename));

		if (outFile.is_open()) {
			outFile << 0 << "\n";
			//Assuming V is nxn
			unsigned int n = sqrt(V.size());
			for (unsigned int i = 0; i < n; i++)
			{
				for (unsigned int j = 0; j < n; j++) {
					outFile << V(i, j) << "\n";
				}
			}
			outFile << 0 << "\n";
		}

		outFile.close();
		return 1;
	}

	void writeToFile(const arma::vec& v, bool overwrite = true) {
		fileOpener(exists(filename));

		if (outFile.is_open()) {
			outFile << 0 << "\n";
			unsigned int n = v.size();
			for (unsigned int i = 0; i < n; i++)
			{
				outFile << v(i) << "\n";
			}
			outFile << 0 << "\n";
			outFile.close();
		}
	}

	inline bool exists(const std::string& name) {
		std::ifstream f(name.c_str());
		return f.good();
	}

	void fileOpener(const bool file) {
		(file) ? outFile.open(filename, std::ofstream::app) : outFile.open(filename, std::ofstream::out);
	}
};