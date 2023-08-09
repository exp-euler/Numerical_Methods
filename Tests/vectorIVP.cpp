#include "TimeIntegration.hpp"
#include <cmath>
#include <getopt.h>

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

d_vector RHS(double t, d_vector y) {
    d_vector rhs(y.size());
    rhs(0) = y(0) + 2*y(1);
    rhs(1) = 3*y(0) + 2*y(1);
    return rhs;
}

d_vector RHS_Lambert1(double t, d_vector y) {
    d_vector rhs(y.size());
    rhs(0) = (-2)*y(0) + y(1) + 2*std::sin(t);
    rhs(1) = y(0) + (-2)*y(1) + 2*(std::cos(t) - std::sin(t));
    return rhs;
}

d_vector RHS_Lambert2(double t, d_vector y) {
    d_vector rhs(y.size());
    rhs(0) = (-2)*y(0) + y(1) + 2*std::sin(t);
    rhs(1) = 998*y(0) + (-999)*y(1) + 999*(std::cos(t) - std::sin(t));
    return rhs;
}

d_vector RHS_Lambert2_nonL(double t, d_vector y) {
    d_vector rhs(y.size());
    rhs(0) = 2*std::sin(t);
    rhs(1) = 999*(std::cos(t) - std::sin(t));
    return rhs;
}

d_vector exact(double t) {
    d_vector sol(2);
    sol(0) = (-8.0/5.0)*std::exp(-t)*(-1) - (4.0/5.0)*std::exp(4*t)*2;
    sol(1) = (-8.0/5.0)*std::exp(-t)*1 - (4.0/5.0)*std::exp(4*t)*3;
    return sol;
}

d_vector exact_Lambert(double t) {
    d_vector sol(2);
    sol(0) = 2*std::exp(-t) + std::sin(t);
    sol(1) = 2*std::exp(-t) + std::cos(t);
    return sol;
}

int main(int argc, char **argv) {

    int steps;

    // Parse for cmd arguments
    char opts[] = ":N:";
    int c = getopt(argc, argv, opts);

    if(c == -1)
        c = 'h';
    do {
        switch (c) {
            case 'N':
                char *endptr;
                steps = std::strtol(optarg, &endptr, 10);

                if(*endptr != '\0') {
                    std::cout << "Steps cannot contain characers other than numbers." << "\n";
                    return 0;
                }

                break;
            case 'h':
                std::cout << "For now, you have to specify number of time steps" << "\n";
                return 0;
            case '?':
                std::cout << "Unrecognized option specified." << "\n";
            case ':':
                std::cout << "Please specify number of time steps." << "\n";
                return 0;
        }
    }while((c = getopt(argc, argv, opts))  != -1);

    // Solve the problem.
    double t0 = 0;
    double t1 = 1;
    double h=(t1-t0)/steps;
    d_vector y0(2);
    y0(0) = 0; y0(1) = -4;

    //TimeIntegration<d_vector> system;
    //system.Solve(ClassicalRK::Euler(), &RHS, h, y0, t0, t1, steps);

    d_vector y0Lam(2);
    y0Lam(0) = 2; y0Lam(1) = 3;

    /*
    TimeIntegration<d_vector> Lambert;
    Lambert.Solve(ClassicalRK::Euler(), &RHS_Lambert2, h, y0Lam, t0, t1, steps);


    Lambert.save_simulation("vector_output.csv", steps);
    */

    d_matrix L(2,2);
    L(0,0)=-2; L(0,1)=1;
    L(1,0)=998; L(1,1)=-999;

    TimeIntegration<d_vector> LambertExp;
    LambertExp.Solve(ExponentialRK<d_matrix>::EEuler(L*(h)), &RHS_Lambert2_nonL, h, y0Lam, t0, t1, steps);

    LambertExp.save_simulation("vector_output.csv", steps);
    return 0;
}
