#ifndef TIMEINTEGRATION
#define TIMEINTEGRATION

#include "Tableaus.hpp"
#include <vector>
#include <functional>

// Provides a function to solve the ODE and stores the approximated solution in Y
// as well as the timesteps T.
template<typename Y_TYPE>
class TimeIntegration
{
private:
    std::vector<Y_TYPE> Y;
    std::vector<double> T;
public:
    void Solve(const ClassicalRK &tableau,
                 std::function<double(double, double)>F, double h, double y0,
                 double t0, double t1, int N);
    // Return by constant reference to avoid both copying and editing.
    const std::vector<double> &getY();
    const std::vector<double> &getT();
};

// ###################################################################
// ###################################################################
//                            IMPLEMENTATIONS
// ###################################################################
// ###################################################################

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
template<typename Y_TYPE>
void TimeIntegration<Y_TYPE>::Solve(const ClassicalRK &tableau,
             std::function<double(double, double)>F, double h, double y0,
             double t0, double t1, int N)
{
    Y.reserve(N);
    T.reserve(N);

    double y = y0;
    double t = t0;

    for(int i=0; i<N; i++)
    {
        std::vector<double> k{F(t,y)};
        for(std::size_t j=0; j<tableau.a.size(); j++)
        {
            double temp = 0;
            for(std::size_t m=0; m<tableau.a[j].size(); m++)
            {
                temp += tableau.a[j][m]*k[m];
            }
            k.push_back(F(t+tableau.c[j]*h, y + temp*h));
        }

        for(std::size_t j=0; j<tableau.b.size(); j++)
        {
            y += tableau.b[j]*h*k[j];
        }
        t += h;

        Y.push_back(y);
        T.push_back(t);
    }

}

template<typename Y_TYPE>
const std::vector<double> &TimeIntegration<Y_TYPE>::getY()
{
    return Y;
}

template<typename Y_TYPE>
const std::vector<double> &TimeIntegration<Y_TYPE>::getT()
{
    return T;
}

#endif
