ODEsolver.o : ODEsolver.cpp ODEsolver.h
	g++ -g -c ODEsolver.cpp

SIRS_test : SIRS_test.cpp ODEsolver.o 
	g++ -g -o SIRS_test SIRS_test.cpp ODEsolver.o 
	
SIRS_vital_dynamics_test : SIRS_vital_dynamics_test.cpp ODEsolver.o 
	g++ -g -o SIRS_vital_dynamics_test SIRS_vital_dynamics_test.cpp ODEsolver.o 
	
SIRS_harmonical_test : SIRS_harmonical_test.cpp ODEsolver.o 
	g++ -g -o SIRS_harmonical_test SIRS_harmonical_test.cpp ODEsolver.o 
	
SIRS_vaccine_test : SIRS_vaccine_test.cpp ODEsolver.o 
	g++ -g -o SIRS_vaccine_test SIRS_vaccine_test.cpp ODEsolver.o 

