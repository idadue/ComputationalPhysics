#include "TridiagonalMatrixSolver.h"

int main()
{
	ThomasSolver t(10, -1, 2, -1);
	t.solve();
	t.display();
	//std::cout << "Time elapsed: " << t.getExecutionTime() << " seconds \n \n";
	//t.writeToFile("test1.csv");

	SpecialThomasSolver sT(10);
	sT.solve();
	sT.display();

	LUSolver LU(10, -1, 2, -1);
	LU.solve();
	LU.display();

	JacobiSolver j(10, 1.0);
	j.analyticalEigenvalues();
	for (int i = 0; i < 10; i++) {
		std::cout << j.getanalyticalEigenvalues()[i] << "\n";
	}

	return 0;
}
