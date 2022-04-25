#ifndef TIMEINTEGRATION
#define TIMEINTEGRATION

#include "Tableaus.hpp"
#include <vector>

// Provides a function to solve the ODE and stores the approximated solution in Y
// as well as the timesteps T.
class TimeIntegration
{
private:
    std::vector<double> Y;
    std::vector<double> T;
public:
    void Solve(const ClassicalRK &tableau,
                 double (*F)(double, double), double h, double y0, double t0,
                 double t1, int N);
    // Return by constant reference to avoid both copying and editing.
    const std::vector<double> &getY();
    const std::vector<double> &getT();
};

#endif
