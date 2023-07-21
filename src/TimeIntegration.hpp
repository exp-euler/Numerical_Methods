//TODO: Pass F by a lambda instead of std::function
//TODO: Add time adjusting capabilities.
#ifndef TIMEINTEGRATION
#define TIMEINTEGRATION

#include "Tableaus.hpp"
#include <fstream>
#include <vector>
#include <functional>
#include <sys/stat.h>

// Provides a function to solve the ODE and stores the approximated
// solution in Y as well as the timesteps T.

// Y_TYPE parameter in case we ever need to work with
// types other than double.
template<typename Y_TYPE>
class TimeIntegration
{
private:
    std::vector<Y_TYPE> Y;
    std::vector<double> T;
public:
    void Solve(const ClassicalRK &tableau,
               std::function<Y_TYPE(double, Y_TYPE)>F, double h,
               Y_TYPE y0, double t0, double t1, int N);
    // Return by constant reference to avoid both copying and editing.
    const std::vector<Y_TYPE> &getY();
    const std::vector<double> &getT();

    void save_simulation(std::string file_name);
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
             std::function<Y_TYPE(double, Y_TYPE)>F, double h, Y_TYPE y0,
             double t0, double t1, int N)
{
    Y.reserve(N);
    T.reserve(N);

    Y_TYPE y = y0;
    double t = t0;
    // Initialize here as work-around to be able to work with
    // scalars and vectors at the same time.
    Y_TYPE temp = y0;

    for(int i=0; i<N; i++)
    {
        std::vector<Y_TYPE> k{F(t,y)};
        for(std::size_t j=0; j<tableau.a.size(); j++)
        {
            temp = 0;
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
const std::vector<Y_TYPE> &TimeIntegration<Y_TYPE>::getY()
{
    return Y;
}

template<typename Y_TYPE>
const std::vector<double> &TimeIntegration<Y_TYPE>::getT()
{
    return T;
}

template<typename Y_TYPE>
void TimeIntegration<Y_TYPE>::save_simulation(std::string file_name)
{
    std::ofstream file_o;

    // Check if file exists.
    struct stat buffer;
    if(stat (file_name.c_str(), &buffer) != 0) {
        file_o.open(file_name);

        // Add a header to the csv file.
        file_o << "time,solution\n";
    }

    // Write the data for each time-step taken.
    for(size_t i=0; i<Y.size(); i++) {
        file_o << T[i] << "," << Y[i] << "\n";
    }

    file_o.close();
}

#endif
