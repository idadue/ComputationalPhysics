#include "sir.h"

SIR::SIR(double h, double N, double S0, double I0, double R0, double a, double b, double c) : h(h), N(N), S0(S0), I0(I0), R0(R0), a(a), b(b), c(c)
{
}

double SIR::dSdt(double S, double I)
{
    return c * (N - S - I) - ((a * S * I) / N) - d * S + e * N - f;
}

double SIR::dIdt(double S, double I)
{
    return ((a * S * I) / N) - b * I - d * I - d_i * I;
}

double SIR::dRdt(double I, double R)
{
    return b * I - c * R - d * R + f;
}

void SIR::output(const std::string &filename, const std::vector<double> &S, const std::vector<double> &I, const std::vector<double> &R)
{
    std::ofstream ofile;
    ofile.open(filename);

    for (int i = 0; i < steps; i++)
    {
        ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
        ofile << std::setw(15) << std::setprecision(8) << S[i];
        ofile << std::setw(15) << std::setprecision(8) << I[i];
        ofile << std::setw(15) << std::setprecision(8) << R[i] << std::endl;
    }
}

double SIR::seasonalVariation(double t, double A, double w, double a_0)
{
    return A * cos(w * t) + a_0;
}

/*

ODESolver Methods below:

*/

void ODESolver::RungeKutta4(double T, const std::string &filename, bool seasonal = false)
{
    steps = (int)(T / h);
    std::cout << T << " " << h << " " << T / h << " " << steps << "\n";
    double k1, l1;
    double k2, l2;
    double k3, l3;
    double k4, l4;

    S.resize(steps, 0);
    I.resize(steps, 0);
    R.resize(steps, 0);

    S[0] = S0;
    I[0] = I0;
    R[0] = R0;

    double initial_f = f;

    if (f_t != 0.0)
    {
        f = 0.0;
    }

    std::cout << a << "\n"
              << b << "\n"
              << c << "\n"
              << d << "\n"
              << d_i << "\n"
              << e << "\n"
              << f << "\n"
              << f_t << "\n";

    for (int i = 0; i < steps; i++)
    {
        if (seasonal)
            SIR::a = seasonalVariation(i * h, 1, 2 * M_PI, 4);

        if (h * i >= f_t && f_t != 0.0)
        {
            //SIR::f = initial_f;
            f = std::min(initial_f, initial_f * 0.1 * (i * h - f_t));
        }

        k1 = h * dSdt(S[i], I[i]);
        l1 = h * dIdt(S[i], I[i]);

        k2 = h * dSdt(S[i] + k1 * 0.5, I[i] + l1 * 0.5);
        l2 = h * dIdt(S[i] + k1 * 0.5, I[i] + l1 * 0.5);

        k3 = h * dSdt(S[i] + k2 * 0.5, I[i] + l2 * 0.5);
        l3 = h * dIdt(S[i] + k2 * 0.5, I[i] + l2 * 0.5);

        k4 = h * dSdt(S[i] + k3, I[i] + l3);
        l4 = h * dIdt(S[i] + k3, I[i] + l3);

        S[i + 1] = S[i] + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        I[i + 1] = I[i] + (1.0 / 6.0) * (l1 + 2 * l2 + 2 * l3 + l4);
        R[i + 1] = N - S[i + 1] - I[i + 1];
    }
    output(filename, S, I, R);
}

void ODESolver::setInitial(double S0, double I0, double R0)
{
    S[0] = S0;
    I[0] = I0;
    R[0] = R0;
}

void ODESolver::setParams(double a, double b, double c, double d, double d_i, double e, double f, double f_t = 0)
{
    SIR::a = a;
    SIR::b = b;
    SIR::c = c;
    SIR::d = d;
    SIR::d_i = d_i;
    SIR::e = e;
    SIR::f = f;
    SIR::f_t = f_t;
}

/*

MCSolver Methods below:

*/

void MCSolver::transitionProbability(int i)
{
    // Probabiltiies of advancing in SIRS
    tProb[0 * 6 + 1] = (a * S[i] * I[i] / total) * h;
    tProb[1 * 6 + 2] = (b * I[i]) * h;
    tProb[2 * 6 + 0] = (c * R[i]) * h;

    // Probabilities of death by unrelated causes
    tProb[0 * 6 + 3] = d * S[i] * h;
    tProb[1 * 6 + 3] = d * I[i] * h;
    tProb[2 * 6 + 3] = d * R[i] * h;

    // Probability of death by infection
    tProb[1 * 6 + 4] = d_i * I[i] * h;

    // Probability of birth
    tProb[5 * 6 + 0] = (e * total) * h;

    // Probability of vaccination
    tProb[0 * 6 + 2] = f * S[i] * h;
}

int main()
{
    ODESolver ode(0.001, 400.0, 300.0, 100.0, 0.0, 4.0, 1.0, 0.5);
    //ode.RungeKutta4(10, "A.txt");
    //ode.setParams(4, 2, 0.5, 0, 0, 0, 0);
    //ode.RungeKutta4(10, "B.txt");
    //ode.RungeKutta4(10, "seasonal.txt", true);

    ode.setParams(4.0, 1.0, 0.5, 0.0, 0.0, 0.0, 80, 70);
    ode.RungeKutta4(100.0, "vaccine.txt", false);

    return 0;
}