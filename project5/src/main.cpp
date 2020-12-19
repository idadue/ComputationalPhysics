#include "sir.h"

int main(int argc, char *argv[])
{
    std::string solver = argv[1];
    double T = std::atof(argv[2]);
    double h = std::atof(argv[3]);
    double a = std::atof(argv[4]);
    double b = std::atof(argv[5]);
    double c = std::atof(argv[6]);
    double d = std::atof(argv[7]);
    double d_I = std::atof(argv[8]);
    double e = std::atof(argv[9]);
    double A = std::atof(argv[10]);
    double w = std::atof(argv[11]) * 2 * M_PI;
    double f = std::atof(argv[12]);
    double f_t = std::atof(argv[13]);
    std::string filename = argv[14];

    if (solver == "rk4")
    {
        ODESolver ode(h, 300.0, 100.0, 0.0, a, b, c);
        ode.setParams(a, b, c, d, d_I, e, f, f_t);

        bool seasonal = false;
        if (A != 0)
        {
            seasonal = true;
            ode.setSeasonalParams(A, w);
        }
        ode.RungeKutta4(T, filename, seasonal);
    }
    else
    {
        printf("MC solver not implemented yet \n");
        exit(0);
    }

    return 0;
}