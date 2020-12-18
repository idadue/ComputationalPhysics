#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>

class SIR
{

public:
    SIR(double h, double N, double S0, double I0, double R0, double a, double b, double c);

    void output(const std::string &filename, const std::vector<double> &S, const std::vector<double> &I, const std::vector<double> &R);

protected:
    double N;  //Population
    double S0; //Initial susceptible
    double I0; //Inital Infected
    double R0; //Inital Recovered

    double a;       //Rate of transmisson
    double b;       //Recovery
    double c;       //Immunity loss
    double d = 0;   //Death rate
    double d_i = 0; //Death rate of infected
    double e = 0;   //Birth rate
    double f = 0;   //Rate of vaccination
    double f_t = 0; //Start vaccination

    double T;
    double h;
    int steps = 0;

    double dSdt(double S, double I);
    double dIdt(double S, double I);
    double dRdt(double I, double R);

    double seasonalVariation(double t, double A_0, double w, double a_0);

private:
};

class ODESolver : public SIR
{
public:
    ODESolver(double h, double N, double S0, double I0, double R0, double a, double b, double c) : SIR(h, N, S0, I0, R0, a, b, c)
    {
        //S.resize(steps, 0);
        //I.resize(steps, 0);
        //R.resize(steps, 0);

        //S[0] = S0;
        //I[0] = I0;
        //R[0] = R0;
        double ok = 0;
    };
    void RungeKutta4(double T, const std::string &filename, bool seasonal);
    void setInitial(double S0, double I0, double R0);
    void setParams(double a, double b, double c, double d, double d_i, double e, double f, double f_t);

private:
    std::vector<double> S;
    std::vector<double> I;
    std::vector<double> R;
};

class MCSolver : public SIR
{
public:
    MCSolver(double steps, double N, double S0, double I0, double R0, double a, double b, double c) : SIR(steps, N, S0, I0, R0, a, b, c)
    {
        S.resize(steps, 0);
        I.resize(steps, 0);
        R.resize(steps, 0);

        S[0] = S0;
        I[0] = I0;
        R[0] = R0;

        total = S0 + I0 + R0;
    }

private:
    double total;
    std::vector<double> S;
    std::vector<double> I;
    std::vector<double> R;
    double tProb[36];

    void transitionProbability(int i);
};