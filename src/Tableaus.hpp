#ifndef TABLEAUS
#define TABLEAUS

#include<vector>

// Provides a way to store an RK method by its Butcher Tableau
// Load the c_i (steps), a_ij (weights) and b_i (final weights).
class ClassicalRK
{
    // All mthods are public for ease of working.
    // Always pass a Tableau Object by const reference to avoid editing tableaus.
    public:
    std::vector<double> c;
    std::vector<std::vector<double>> a;
    std::vector<double> b;

    ClassicalRK(std::vector<double> steps,
                std::vector<std::vector<double>> weights,
                std::vector<double> weightsHigher);
};

// Declare (only!) the methods tableaus below before defining them in Tableaus.cpp
extern ClassicalRK Euler;
extern ClassicalRK RK4;

#endif //TABLEAUS
