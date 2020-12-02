// An example of the use of methods in the ODEsolver class to numerically
// solve an SIR model

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
  while (time<=(10.0+delta_t)){
     vtime.push_back(time);
     time = time + delta_t;
  }
  
  // set parameters and initial values
  double a = 4.0;
  double b = 1.0;
  double c = 0.5;
  vector<double> Parameters;
  Parameters.push_back(a);
  Parameters.push_back(b);
  Parameters.push_back(c);

  double N = 400;
  double I_0 = 1.0;
  double S_0 = N-I_0;
  double R_0 = 0.0;
  vector<double> InitialValues;
  InitialValues.push_back(S_0);
  InitialValues.push_back(I_0);
  InitialValues.push_back(R_0);
  
  // ModelType==1 is the SIRS model
  int ModelType = 1;
  ODEsolver SIRS(InitialValues, Parameters, ModelType);

  // solve the model using the 4th order Runge-Kutta method
  vector<vector<double> > ModelEstimates;
  SIRS.SimulateModel(vtime
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

