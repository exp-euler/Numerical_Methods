#include "TimeIntegration.hpp"
#include "Tableaus.hpp"
#include "Vector.hpp"
#include <cstdlib>
#include <vector>
#include <iostream>
#include <time.h>

#include <unistd.h>

// Demonstration of specifying a Runge-Kutta method to solve an equation like:
//
//                     dy/dt = F(t,y),    y(0) = y0
//
// where the function F(t,y) can be arbitrary and is referred as Right Hand Side.
double RHS(double t, double y)
{
    return 1+t;
}

int main(int argc, char*argv[])
{
    // Testing of the time integration part of the code (serial for now)
    TimeIntegration<double> equation;
    equation.Solve(ClassicalRK::Euler(), &RHS, (0.0+1.0)/100, 2.0, 0, 1, 100);
    //const std::vector<double> &sol = equation.getY();

    equation.save_simulation("scalar_output.csv", 100);
    return 0;
}
