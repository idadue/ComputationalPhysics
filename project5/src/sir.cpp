#include "sir.h"

SIR::SIR(double h, double S0, double I0, double R0, double a, double b, double c) : h(h), S0(S0), I0(I0), R0(R0), a(a), b(b), c(c)
{
}

double SIR::dSdt(double S, double I, double R, double N)
{
    return c * R - ((a * S * I) / N) - d * S + e * N - f;
}

double SIR::dIdt(double S, double I, double N)
{
    return ((a * S * I) / N) - b * I - d * I - d_i * I;
}

double SIR::dRdt(double R, double I)
{
    return b * I - c * R - d * R + f;
}

void SIR::output(const std::string &filename, const std::vector<double> &S, const std::vector<double> &I, const std::vector<double> &R, const std::vector<double> &N)
{
    std::ofstream ofile;
    ofile.open(filename);

    for (int i = 0; i < steps; i++)
    {
        ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
        ofile << std::setw(15) << std::setprecision(8) << S[i];
        ofile << std::setw(15) << std::setprecision(8) << I[i];
        ofile << std::setw(15) << std::setprecision(8) << R[i];
        ofile << std::setw(15) << std::setprecision(8) << N[i] << std::endl;
    }
}

double SIR::seasonalVariation(double t, double A, double w, double a_0)
{
    return A * cos(w * t) + a_0;
}

/*

ODESolver Methods below:

*/

/*Implementation of Runge Kutta 4th order method.*/
void ODESolver::RungeKutta4(double T, const std::string &filename, bool seasonal = false)
{
    steps = (int)(T / h);
    std::cout << T << " " << h << " " << steps << "\n";
    double k1, l1, p1;
    double k2, l2, p2;
    double k3, l3, p3;
    double k4, l4, p4;

    S.resize(steps, 0);
    I.resize(steps, 0);
    R.resize(steps, 0);
    N.resize(steps, 0);

    S[0] = S0;
    I[0] = I0;
    R[0] = R0;
    N[0] = S0 + I0 + R0;

    double initial_f = f;

    std::cout << "a = " << a << "\n"
              << "b = " << b << "\n"
              << "c = " << c << "\n"
              << "d = " << d << "\n"
              << "d_i = " << d_i << "\n"
              << "e = " << e << "\n"
              << "f = " << f << "\n"
              << "f_t = " << f_t << "\n";

    if (f_t != 0.0)
    {
        f = 0.0;
    }

    for (int i = 0; i < steps; i++)
    {
        if (seasonal)
            SIR::a = seasonalVariation(i * h, A, w, a_0);

        if (h * i >= f_t && f_t != 0.0)
        {
            f = std::min(initial_f, initial_f * 0.1 * (i * h - f_t));
        }

        k1 = h * dSdt(S[i], I[i], R[i], N[i]);
        l1 = h * dIdt(S[i], I[i], N[i]);
        p1 = h * dRdt(R[i], I[i]);

        k2 = h * dSdt(S[i] + k1 * 0.5, I[i] + l1 * 0.5, R[i] + p1 * 0.5, N[i]);
        l2 = h * dIdt(S[i] + k1 * 0.5, I[i] + l1 * 0.5, N[i]);
        p2 = h * dRdt(R[i] + p1 * 0.5, I[i] + l1 * 0.5);

        k3 = h * dSdt(S[i] + k1 * 0.5, I[i] + l1 * 0.5, R[i] + p2 * 0.5, N[i]);
        l3 = h * dIdt(S[i] + k2 * 0.5, I[i] + l2 * 0.5, N[i]);
        p3 = h * dRdt(R[i] + p2 * 0.5, I[i] + l2 * 0.5);

        k4 = h * dSdt(S[i] + k3, I[i] + l3, R[i] + p3, N[i]);
        l4 = h * dIdt(S[i] + k3, I[i] + l3, N[i]);
        p4 = h * dRdt(R[i] + p3, I[i] + l3);

        S[i + 1] = S[i] + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        I[i + 1] = I[i] + (1.0 / 6.0) * (l1 + 2 * l2 + 2 * l3 + l4);
        R[i + 1] = R[i] + (1.0 / 6.0) * (p1 + 2 * p2 + 2 * p3 + p4);

        N[i + 1] = S[i + 1] + I[i + 1] + R[i + 1];
    }
    output(filename, S, I, R, N);
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

void ODESolver::setSeasonalParams(double A, double w)
{
    SIR::A = A;
    SIR::w = w;
    SIR::a_0 = a;
}

/*

MCSolver Methods below:

*/