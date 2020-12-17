//#######################################################################
// Functions for solving the differential equations using the 4th order 
// Runge-Kutta method
//
// This class have been implemented using the excellent example of a SIR 
// model provided by:
// Author: Sherry Towers (smtowers@asu.edu), Copyright Sherry Towers, 2013
//########################################################################

#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include "ODEsolver.h"    

using namespace std;

//########################################################################
// InitialValues are the initial starting values for the compartments.
// Parameters    is the vector of parameters for the model.  
// ModelType     1 is a standard SIRS model
//               2 is a SIRS model with vital dynamics
//               3 is a SIRS model with harmonic transmission
//               4 is a vaccination model. 
//########################################################################
   ODEsolver::ODEsolver(vector<double> InitialValues
                           ,vector<double> Parameters
                           ,int            ModelType)
                           :_InitialValues(InitialValues)
                           ,_Parameters(Parameters)
                           ,_ModelType(ModelType){
   }

//########################################################################
// For _ModelType==1 (SIRS model)
//    Compartments are (in order) S I R
//    and the parameters in _Parameters are (in order) a b c
// For _ModelType==2 (SIRS model with vital dynamics)
//    Compartments are (in order) S I R
//    and the parameters in _Parameters are (in order) a b c d d_i e
// For _ModelType==3 (harmonic SIRS model)
//    Compartments are (in order) S I R
//    and the parameters in _Parameters are (in order) a0 b c A omega
//  For _ModelType==4 (vaccination model)
//    Compartments are (in order) S I R
//    and the parameters in _Parameters are (in order) a b c f
//########################################################################
   void ODEsolver::ModelEquations(double          time
                                   ,vector<double>  Compartments
                                   ,vector<double>& Derivatives){
     Derivatives.clear();
     
     double N = accumulate(Compartments.begin()
                          ,Compartments.end()
                          ,0.0
                          ,plus<double>());

     //###################################################################
     // SIRS model
     //   Compartments[0] = S
     //   Compartments[1] = I
     //   Compartments[2] = R
     //
     //   _Parameters[0] = a
     //   _Parameters[1] = b
     //   _Parameters[2] = c 
     //###################################################################
     if (_ModelType==1){
       double S = Compartments[0];
       double I = Compartments[1];
       double R = Compartments[2];
       double a  = _Parameters[0];
       double b = _Parameters[1];
       double c = _Parameters[2];

       Derivatives.push_back(c*R-a*S*I/N);
       Derivatives.push_back(a*S*I/N - b*I);
       Derivatives.push_back(b*I-c*R);

     //###################################################################
     // SIRS model with vital dynamics
     //   Compartments[0] = S
     //   Compartments[1] = I
     //   Compartments[2] = R
     //
     //   _Parameters[0] = a
     //   _Parameters[1] = b
     //   _Parameters[2] = c
     //   _Parameters[3] = d
     //   _Parameters[4] = d_i
     //   _Parameters[5] = e
     //###################################################################
     }else if (_ModelType==2){
       double S = Compartments[0];
       double I = Compartments[1];
       double R = Compartments[2];
       double a  = _Parameters[0];
       double b = _Parameters[1];
       double c = _Parameters[2];
       double d = _Parameters[3];
       double d_i = _Parameters[4];
       double e = _Parameters[5];

       Derivatives.push_back(c*R-a*S*I/N-d*S+e*N);
       Derivatives.push_back(a*S*I/N-b*I-d*I-d_i*I);
       Derivatives.push_back(b*I-c*R-d*R);

     //###################################################################
     // Harmonic SIRS to model seasonal variations
     //   Compartments[0] = S
     //   Compartments[1] = I
     //   Compartments[2] = R
     //
     //   _Parameters[0] = a0
     //   _Parameters[1] = b
     //   _Parameters[2] = c
     //   _Parameters[3] = A
     //   _Parameters[4] = omega
     //###################################################################
     }else if (_ModelType==3){
       double S = Compartments[0];
       double I = Compartments[1];
       double R = Compartments[2];
       double a0  = _Parameters[0];
       double b   = _Parameters[1];
       double c   = _Parameters[2];
       double A   = _Parameters[3];
       double omega = _Parameters[4];       
       
       double a = a0*(1.0+A*cos(omega*time));

       Derivatives.push_back(c*R-a*S*I/N);
       Derivatives.push_back(a*S*I/N - b*I);
       Derivatives.push_back(b*I-c*R);

     //###################################################################
     // Vaccination
     //   Compartments[0] = S
     //   Compartments[1] = I
     //   Compartments[2] = R
     //
     //   _Parameters[0] = a
     //   _Parameters[1] = b
     //   _Parameters[2] = c
     //   _Parameters[3] = f
     //###################################################################
     }else if (_ModelType==4){
       double S = Compartments[0];
       double I = Compartments[1];
       double R = Compartments[2];
       double a  = _Parameters[0];
       double b = _Parameters[1];
       double c = _Parameters[2];
       double f = _Parameters[3];
       double vaccineTime = _Parameters[4];
       
       Derivatives.push_back(c*R-a*S*I/N-f);
       Derivatives.push_back(a*S*I/N-b*I);
       Derivatives.push_back(b*I-c*R+f);
     }

     if (Derivatives.size()!=Compartments.size()){
       cout << "ODEsolver(Compartments): Number of derivatives must match number of compartments!\n";
       cout << "                       Number of derivatives  = " << Derivatives.size() << endl;
       cout << "                       Number of compartments = " << Compartments.size() << endl;
     }
   }

   void ODEsolver::SimulateModel(vector<double>  vtime
                                  ,vector<vector<double> >& ModelEstimates
                                  ){
     vector<vector<double> > ModelEstimatesTranspose;
     SimulateModel(vtime
                  ,ModelEstimates
                  ,ModelEstimatesTranspose);
   }


   void ODEsolver::SimulateModel(vector<double>  vtime
                                  ,vector<vector<double> >& ModelEstimates
                                  ,vector<vector<double> >& ModelEstimatesTranspose
                                  ){
     ModelEstimates.clear();
     ModelEstimatesTranspose.clear();

     double newI;
     vector<double> atemp = _InitialValues;
     vector<double> oldModelEstimates = _InitialValues;
     if (_ModelType==1){ // tack on the number of new infections per unit time a*S*I/N
       double N = accumulate(atemp.begin()
                            ,atemp.end()
                            ,0.0
                            ,plus<double>());

       newI = _Parameters[0]*atemp[0]*atemp[1]/N;
     }
     if (_ModelType==2){ // tack on the number of new infections per unit time kappa*E
       newI = _Parameters[2]*atemp[1];
     }

     // set up the ModelEstimates and ModelEstimatesTranspose vectors
     if (_ModelType<=2) atemp.push_back(newI);
     ModelEstimates.push_back(atemp);

     for (int i=0;i<atemp.size();i++){
       vector<double> btemp;
       btemp.push_back(atemp[i]);
       ModelEstimatesTranspose.push_back(btemp);
     }

     // simulation over the time steps
     for (int i=1;i<vtime.size();i++){
        double delta_t = vtime[i]-vtime[i-1];
        vector<double> newModelEstimates;
        RungeKutta4(delta_t,vtime[i],oldModelEstimates,newModelEstimates);
        oldModelEstimates = newModelEstimates;

        vector<double> newModelEstimatesb = newModelEstimates;
        if (_ModelType==1 || _ModelType==2 || _ModelType==4 ){ // tack on the number of new infections per unit time a*S*I/N
          double N = accumulate(newModelEstimates.begin()
                               ,newModelEstimates.end()
                               ,0.0
                               ,plus<double>());

          double newI = _Parameters[0]*newModelEstimates[0]*newModelEstimates[1]/N;
          newModelEstimatesb.push_back(newI);
        }
       
        ModelEstimates.push_back(newModelEstimatesb);

        for (int j=0;j<newModelEstimatesb.size();j++){
          ModelEstimatesTranspose[j].push_back(newModelEstimatesb[j]);
        }
     }
   }

   void ODEsolver::RungeKutta4(double  delta_t   
                                     ,double  time
                                     ,vector<double>  oldModelEstimates
                                     ,vector<double>& newModelEstimates){
      newModelEstimates.clear();
      vector<double> updatedModelEstimates;
      vector<double> Derivatives;
      vector<double> a1;
      ModelEquations(time,oldModelEstimates,Derivatives);
      for (int i=0;i<oldModelEstimates.size();i++){
         a1.push_back(delta_t*Derivatives[i]); 
         updatedModelEstimates.push_back(oldModelEstimates[i]+delta_t*Derivatives[i]/2.0); 
      }

      ModelEquations(time,updatedModelEstimates,Derivatives);
      updatedModelEstimates.clear();
      vector<double> a2;
      for (int i=0;i<oldModelEstimates.size();i++){
         a2.push_back(delta_t*Derivatives[i]); 
         updatedModelEstimates.push_back(oldModelEstimates[i]+delta_t*Derivatives[i]/2.0); 
      }

      ModelEquations(time,updatedModelEstimates,Derivatives);
      updatedModelEstimates.clear();
      vector<double> a3;
      for (int i=0;i<oldModelEstimates.size();i++){
         a3.push_back(delta_t*Derivatives[i]); 
         updatedModelEstimates.push_back(oldModelEstimates[i]+delta_t*Derivatives[i]); 
      }

      vector<double> a4;
      ModelEquations(time,updatedModelEstimates,Derivatives);
      for (int i=0;i<oldModelEstimates.size();i++){
         a4.push_back(delta_t*Derivatives[i]); 
      }

      for (int i=0;i<oldModelEstimates.size();i++){
        newModelEstimates.push_back(oldModelEstimates[i]
                                   +(a1[i]
                                    +2.0*a2[i]
                                    +2.0*a3[i]
                                    +a4[i])/6.0);
        if (newModelEstimates[i]<0) newModelEstimates[i] = 0;
      }

   }

   ODEsolver::~ODEsolver(){
   }

