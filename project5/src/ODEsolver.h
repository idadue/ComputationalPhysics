//#######################################################################
// Functions for solving the differential equations using the 4th order 
// Runge-Kutta method
//
// This class have been implemented using the excellent example of a SIR 
// model provided by:
// Author: Sherry Towers (smtowers@asu.edu), Copyright Sherry Towers, 2013
//########################################################################
#include <iostream>
#include <vector>

class ODEsolver{
private:
  std::vector<double> _InitialValues;
  std::vector<double> _Parameters;
  int                 _ModelType;

//########################################################################
//  ModelType==1 is an SIRS model
//    _IntialValues in order S, I, R
//    _Parameters   in order a, b, c
//
//  ModelType==2 is an SIRS model with vital dynamics
//    _IntialValues in order S, I, R
//    _Parameters   in order a, b, c, d, d_i, e
//
//  ModelType==3 is a harmonic SIRS model
//    _IntialValues in order S, I, R
//    _Parameters   in order a0, b, c, A, omega
//
//  ModelType==4 is a vaccination model
//    _IntialValues in order S, I, R
//    _Parameters   in order a, b, c, f
//########################################################################
   void ModelEquations(double time, std::vector<double>  ModelEstimates
                      ,std::vector<double>& Derivatives);

// 4 th order Runge-Kutta mehtod
   void RungeKutta4(double time_step, double time, std::vector<double> oldModelEstimates
                        ,std::vector<double>& newModelEstimates);

public:
   ODEsolver(std::vector<double> InitialValues
          ,std::vector<double> Parameters
          ,int                 ModelType);
 
   ~ODEsolver();
 
   void SimulateModel(std::vector<double>  vtime
                     ,std::vector<std::vector<double> >& ModelEstimates
                     ,std::vector<std::vector<double> >& ModelEstimatesTranspose
                     );

   void SimulateModel(std::vector<double>  vtime
                     ,std::vector<std::vector<double> >& ModelEstimates
                     );
};

