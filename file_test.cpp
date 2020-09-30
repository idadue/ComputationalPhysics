#include <iostream>
#include <fstream>


int main() {

	std::fstream ofile;
	ofile.open("...//Results//haha.csv", std::fstream::out);

	ofile << "Works?" << std::endl;

	ofile.close();
	return 0;
}