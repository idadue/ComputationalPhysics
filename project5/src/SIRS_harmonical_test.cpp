// An example of the use of methods in the ODEsolver class to numerically
// solve an SIRS model

#include <algorithm>
#include <cstdlib>
#include <vector>
#include <iostream>
#include "ODEsolver.h"

using namespace std;
int main(int argc, char **argv){

  double time = 0.0;
  double delta_t = 0.0010;
  vector<double> vtime;
  while (time<=(100.0+delta_t)){
     vtime.push_back(time);
     time = time + delta_t;
  }
  
  // set parameters and initial values
  double a0  = 4.0;
  double b = std::atoi(argv[1]);
  double c = 0.5;
  double A = std::atoi(argv[2]);
  double omega = 0.05*2.0*3.14;
  
  vector<double> Parameters;
  Parameters.push_back(a0);
  Parameters.push_back(b);
  Parameters.push_back(c);
  Parameters.push_back(A);
  Parameters.push_back(omega);

  double N = 400;
  double I_0 = 100.0;
  double S_0 = N-I_0;
  double R_0 = 0.0;
  vector<double> InitialValues;
  InitialValues.push_back(S_0);
  InitialValues.push_back(I_0);
  InitialValues.push_back(R_0);
  
  // ModelType==3 is the SIRS model with seasonal variations
  int ModelType = 3;
  ODEsolver SIRS_harmonical(InitialValues, Parameters, ModelType);

  // solve the model using the 4th order Runge-Kutta method
  vector<vector<double> > ModelEstimates;
  SIRS_harmonical.SimulateModel(vtime
                     ,ModelEstimates);

  for (int i=0;i<vtime.size();i++){
     cout << vtime[i] << " ";
     for (int j=0;j<InitialValues.size();j++){
        cout << ModelEstimates[i][j] << " ";
     }
     cout << endl;
  }
 
  return 0;

}

