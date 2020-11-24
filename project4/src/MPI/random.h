#include <random>

class Random
{

public:
    Random();
    Random(double seed);
    double generateMersenneTwisterNumber(int min, int max);
    int generateMersenneTwistNumber(int min, int max);
    void setSeed(double seed);
    int spin();

private:
    double seed;
    std::mt19937 gen;
};