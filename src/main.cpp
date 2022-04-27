#include "TimeIntegration.hpp"
#include "Tableaus.hpp"
#include <vector>

// Demonstration of specifying a Runge-Kutta method to solve an equation like:
//
//                     dy/dt = F(t,y),    y(0) = y0
//
// where the function F(t,y) can be arbitrary and is referred as Right Hand Side.
double RHS(double t, double y)
{
    return 1+t;
}

int main()
{
    TimeIntegration equation;
    equation.Solve(ClassicalRK::Euler(), &RHS, (0.0+1.0)/100, 2.0, 0, 1, 100);
    const std::vector<double> &sol = equation.getY();

    return 0;
}
