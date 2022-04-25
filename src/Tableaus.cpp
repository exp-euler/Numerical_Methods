#include "Tableaus.hpp"
#include <vector>

ClassicalRK::ClassicalRK(std::vector<double> steps,
             std::vector<std::vector<double>> weights,
             std::vector<double> weightsHigher)
{
    c = steps;
    a = weights;
    b = weightsHigher;
}


// Define the method tableaus below (after declaring in Tableaus.cpp !)
ClassicalRK Euler({},{},{1});
ClassicalRK RK4({1.0/2, 1.0/2, 1},
                {{1.0/2}, {0, 1.0/2}, {0, 0, 1}},
                {1.0/6, 1.0/3, 1.0/3, 1.0/6});
