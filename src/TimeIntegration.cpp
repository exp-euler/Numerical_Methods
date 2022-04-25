#include "TimeIntegration.hpp"
#include <vector>
#include <cmath>

// Solver function that allows the user to specify:
//      the method by              tableau
//      the right hand side by     F
//      step size by               h
//      initial value by           y0
//      initial time by            t0
//      final time by              t1
//      number of steps by         N
//
// For the algorithm used inside the for loop, please see the Wiki page
// specified in the README file.
void TimeIntegration::Solve(const ClassicalRK &tableau,
             double (*F)(double, double), double h, double y0, double t0,
             double t1, int N)
{
    double y = y0;
    double t = t0;

    for(int i=0; i<N; i++)
    {
        std::vector<double> k{F(t,y)};
        for(unsigned int j=0; j<tableau.a.size(); j++)
        {
            double temp = 0;
            for(unsigned int m=0; m<tableau.a[j].size(); m++)
            {
                temp += tableau.a[j][m]*k[m];
            }
            k.push_back(F(t+tableau.c[j]*h, y + temp*h));
        }

        for(unsigned int j=0; j<tableau.b.size(); j++)
        {
            y += tableau.b[j]*h*k[j];
        }
        t += h;

        Y.push_back(y);
        T.push_back(t);
    }

}

const std::vector<double> &TimeIntegration::getY()
{
    return Y;
}

const std::vector<double> &TimeIntegration::getT()
{
    return T;
}
