//TODO: Pass F by a lambda instead of std::function
//TODO: Add time adjusting capabilities.
#ifndef TIMEINTEGRATION
#define TIMEINTEGRATION

#include "Matrix.hpp"
#include "Tableaus.hpp"
#include <fstream>
#include <iomanip>
#include <vector>
#include <functional>
#include <sys/stat.h>

#ifdef EIGEN_YES

#include <Eigen/Dense>
typedef Eigen::VectorXd d_vector;
typedef Eigen::MatrixXd d_matrix;

#else

#include "Matrix.hpp"
#include "Vector.hpp"
typedef Vector<double> d_vector;
typedef Matrix<double> d_matrix;

#endif

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

    // TODO: Find a better way to do this!
    double mult(double &a, double &b);
    d_vector mult(d_vector &D, d_vector &V);
    d_vector mult(d_matrix &M, d_vector &V);

    void file_stream(std::ofstream &file_o, d_vector &placeholder);
    void file_stream(std::ofstream &file_o, double placeholder);
public:
    template<typename TABLEAU>
    void Solve(TABLEAU tableau,
               std::function<Y_TYPE(double, Y_TYPE)>F, double h,
               Y_TYPE y0, double t0, double t1, int N);
    // Return by constant reference to avoid both copying and editing.
    const std::vector<Y_TYPE> &getY();
    const std::vector<double> &getT();

    void save_simulation(std::string file_name, int tau);
};

// ###################################################################
// ###################################################################
//                            IMPLEMENTATIONS
// ###################################################################
// ###################################################################

template<typename Y_TYPE>
double TimeIntegration<Y_TYPE>::mult(double &a, double &b) {
    return a*b;
}

template<typename Y_TYPE>
d_vector TimeIntegration<Y_TYPE>::mult(d_vector &D, d_vector &V) {
    return D.cwiseProduct(V);
}

template<typename Y_TYPE>
d_vector TimeIntegration<Y_TYPE>::mult(d_matrix &M, d_vector &V) {
    return M*V;
}

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
template<typename TABLEAU>
void TimeIntegration<Y_TYPE>::Solve(TABLEAU tableau,
             std::function<Y_TYPE(double, Y_TYPE)>F, double h, Y_TYPE y0,
             double t0, double t1, int N)
{
    Y.reserve(N);
    T.reserve(N);

    Y_TYPE y = y0;
    double t = t0;
    // Initialize here as work-around to be able to work with
    // scalars and vectors at the same time.
    Y_TYPE intermediate_y = y0;

    for(int i=0; i<N; i++)
    {
        std::vector<Y_TYPE> k{F(t,y)};
        for(std::size_t j=0; j<tableau.a.size(); j++)
        {
            intermediate_y = mult(tableau.e[j], y);
            for(std::size_t m=0; m<tableau.a[j].size(); m++)
            {
                intermediate_y += h*mult(tableau.a[j][m],k[m]);
            }
            k.push_back(F(t+ h*tableau.c[j],  intermediate_y));
        }

        y = mult(tableau.e[tableau.b.size()-1],y);
        for(std::size_t j=0; j<tableau.b.size(); j++)
        {
            y = y + h*mult(tableau.b[j],k[j]);
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

// TODO: Find a better way to achieve writing to file for scalar and vector problems
template<typename Y_TYPE>
void TimeIntegration<Y_TYPE>::file_stream(std::ofstream &file_o, d_vector &placeholder){
    // Write the data for each time-step taken.
    for(size_t i=0; i<Y.size(); i++) {
        file_o << std::setprecision(std::numeric_limits<double>::max_digits10) << T[i] << "," ;
        for(size_t j=0; j<Y[0].size(); j++)
            file_o << std::setprecision(std::numeric_limits<double>::max_digits10) << Y[i](j) << "," ;

        file_o << std::setprecision(std::numeric_limits<double>::max_digits10) << "\n" ;
    }

}

// TODO: Find a better way to achieve writing to file for scalar and vector problems
template<typename Y_TYPE>
void TimeIntegration<Y_TYPE>::file_stream(std::ofstream &file_o, double placeholder){
    // Write the data for each time-step taken.
    for(size_t i=0; i<Y.size(); i++) {
        file_o << std::setprecision(std::numeric_limits<double>::max_digits10) << T[i] << "," << Y[i] << "," << "\n";
    }

}

template<typename Y_TYPE>
void TimeIntegration<Y_TYPE>::save_simulation(std::string file_name, int tau)
{
    file_name.erase(file_name.end()-4, file_name.end());
    file_name = file_name + "_" + std::to_string(tau) + ".csv";
    std::ofstream file_o;

    // Check if file exists.
    struct stat buffer;
    if(stat (file_name.c_str(), &buffer) != 0) {
        file_o.open(file_name);

        // Add a header to the csv file.
        //file_o << "time,";
        //for(size_t i=0; i<2; i++)
        //    file_o << "x"<< i << ",";
        //file_o << "\n";
    }

    TimeIntegration::file_stream(file_o, Y[0]);

    file_o.close();
}

#endif
