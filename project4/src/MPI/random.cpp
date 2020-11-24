#include "random.h"

Random::Random()
{
    std::random_device rd;
    seed = rd();
    gen.seed(seed);
}

Random::Random(double seed) : seed(seed)
{
    gen.seed(seed);
}

double Random::generateMersenneTwisterNumber(int min, int max)
{
    std::uniform_real_distribution<double> distrubution(min, max);
    return distrubution(gen);
}

int Random::generateMersenneTwistNumber(int min, int max)
{
    std::uniform_int_distribution<int> distrubution(min, max);
    return distrubution(gen);
}

void Random::setSeed(double seed)
{
    this->seed = seed;
}

int Random::spin()
{
    /*
    Return a either a -1 or a 1. Used for Ising model simulations.
    */
    double num = generateMersenneTwisterNumber(0, 1);
    return (num > 0.5) ? 1 : -1;
}
