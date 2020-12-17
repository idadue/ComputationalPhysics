// An example of the use of methods in the ODEsolver class to numerically
// solve an SIRS model with vital dynamics

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
  while (time<=(40.0+delta_t)){
     vtime.push_back(time);
     time = time + delta_t;
  }
  
  // set parameters and initial values
  double a  = 4.0;
  double b = std::atoi(argv[1]);
  double c = 0.5;
  double d = 0.003;
  double d_i = 0.01;
  double e = 0.0005;
  vector<double> Parameters;
  Parameters.push_back(a);
  Parameters.push_back(b);
  Parameters.push_back(c);
  Parameters.push_back(d);
  Parameters.push_back(d_i);
  Parameters.push_back(e);

  double N = 400;
  double I_0 = 100.0;
  double S_0 = N-I_0;
  double R_0 = 0.0;
  vector<double> InitialValues;
  InitialValues.push_back(S_0);
  InitialValues.push_back(I_0);
  InitialValues.push_back(R_0);
  
  // ModelType==2 is the SIRS model with vital dynamics
  int ModelType = 2;
  ODEsolver SIRS_vital(InitialValues, Parameters, ModelType);

  // solve the model using the 4th order Runge-Kutta method
  vector<vector<double> > ModelEstimates;
  SIRS_vital.SimulateModel(vtime
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

