#include <random>

class Random
{
    /*
    Simple class for finding random numbers.
    Currently only supports number generated by the 
    Mersienne-Twister engine. 
    */

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