#include "disease.h"
#include <math.h>

int main(int argc, char *argv[])
{
    double t_total = std::atoi(argv[1]);
    double a = std::atof(argv[2]);
    double b = std::atof(argv[3]);
    double c = std::atof(argv[4]);
    double d = std::atof(argv[5]);
    double d_I = std::atof(argv[6]);
    double e = std::atof(argv[7]);
    double A = std::atof(argv[8]);
    double w = std::atof(argv[9]) * 2 * M_PI;
    double f = std::atof(argv[10]);
    double fT = std::atof(argv[11]);
    std::string filename = argv[12];

    disease dis;
    dis.initialize(t_total, 300, 100, 0, a, b, c, d, d_I, e, A, w, f, fT);

    dis.solver(filename);

    return 0;
}