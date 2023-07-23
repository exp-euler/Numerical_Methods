#include "Tableaus.hpp"
#include <vector>

ClassicalRK::ClassicalRK(std::vector<double> steps,
             std::vector<double> t_factors,
             std::vector<std::vector<double>> weights,
             std::vector<double> weightsHigher)
{
    c.reserve(steps.size());
    e.reserve(t_factors.size());
    a.reserve(weights.size());
    b.reserve(weightsHigher.size());

    c = steps;
    e = t_factors;
    a = weights;
    b = weightsHigher;
}

ClassicalRK ClassicalRK::Euler()
{
    return ClassicalRK({},{1},{},{1});
}

ClassicalRK ClassicalRK::RK4()
{
    return ClassicalRK({1.0/2, 1.0/2, 1},
                       {1, 1, 1, 1},
                       {{1.0/2}, {0, 1.0/2}, {0, 0, 1}},
                       {1.0/6, 1.0/3, 1.0/3, 1.0/6});
}
