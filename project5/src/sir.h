#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>

class SIR
{

public:
    SIR(double h, double S0, double I0, double R0, double a, double b, double c);

    void output(const std::string &filename, const std::vector<double> &S, const std::vector<double> &I, const std::vector<double> &R, const std::vector<double> &N);

protected:
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

    double A = 0;
    double w = 0;
    double a_0;

    double T;
    double h;
    int steps = 0;

    double dSdt(double S, double I, double R, double N);
    double dIdt(double S, double I, double N);
    double dRdt(double R, double I);

    double seasonalVariation(double t, double A_0, double w, double a_0);

private:
};

class ODESolver : public SIR
{
public:
    ODESolver(double h, double S0, double I0, double R0, double a, double b, double c) : SIR(h, S0, I0, R0, a, b, c){};
    void RungeKutta4(double T, const std::string &filename, bool seasonal);
    void setInitial(double S0, double I0, double R0);
    void setParams(double a, double b, double c, double d, double d_i, double e, double f, double f_t);
    void setSeasonalParams(double A, double w);

private:
    std::vector<double> S;
    std::vector<double> I;
    std::vector<double> R;
    std::vector<double> N;
};

class MCSolver : public SIR
{
public:
    MCSolver(double steps, double N, double S0, double I0, double R0, double a, double b, double c) : SIR(steps, S0, I0, R0, a, b, c)
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