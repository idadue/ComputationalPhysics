#include "TridiagonalMatrixSolver.h"

int main()
{
	ThomasSolver t(100000, -1, 2, -1);
	t.solve();
	//t.display();
	std::cout << "Time elapsed: " << t.getExecutionTime() << " seconds \n \n";
	//t.writeToFile("test1.csv");

	SpecialThomasSolver sT(100);
	sT.solve();
	sT.display();

	LUSolver LU(100, -1, 2, -1);
	LU.solve();
	LU.display();

	return 0;
}