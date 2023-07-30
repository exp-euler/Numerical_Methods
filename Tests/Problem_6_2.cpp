#include "Matrix.hpp"
#include "Tableaus.hpp"
#include "TimeIntegration.hpp"
#include "Vector.hpp"
#include <cmath>
#include <getopt.h>

typedef Vector<double> d_vector;
typedef Matrix<double> d_matrix;

// TODO: diffusion coefficient currently taken to be K=1
d_vector exact_P62_sin(double t, int Psize) {
    d_vector sol(Psize);

    double x;
    double dx = (1.0-0.0)/(Psize+1);
    for(int i=0; i<Psize; i++) {
        x = dx * (i+1);
        sol[i] = 10*(1 - x)*x*(1+std::sin(t)) + 2;
    }

    return sol;
}

d_vector RHS_P62_sin(double t, d_vector y) {
    d_vector rhs(y.Size());

    double x;
    double dx = (1.0-0.0)/(y.Size()+1);
    for(int i=0; i<y.Size(); i++) {
        x = dx * (i+1);
        rhs[i] = 10*(1 - x)*x*(std::cos(t)) + 2*10*(1+std::sin(t)) - 1/(1+(x*(1-x)*(10+10*std::sin(t))+2) * (x*(1-x)*(10+10*std::sin(t))+2)) // Phi(x,t)
        + 1/(1+y[i]*y[i]); // Fraction 1/(1+y^2)
    }

    // Boundary condition modification specific to Laplacian 1 -2 1
    rhs[0] += (1/(dx*dx))*2; rhs[y.Size()-1] += (1/(dx*dx))*2;

    return rhs;
}

d_vector RHS_P62_sin_wL(double t, d_vector y) {
    d_vector rhs(y.Size());

    double x;
    double dx = (1.0-0.0)/(y.Size()+1);
    for(int i=0; i<y.Size(); i++) {
        x = dx * (i+1);
        rhs[i] = 10*(1 - x)*x*(std::cos(t)) + 2*10*(1+std::sin(t)) - 1/(1+(x*(1-x)*(10+10*std::sin(t))+2) * (x*(1-x)*(10+10*std::sin(t))+2))
        + 1/(1+y[i]*y[i]);
    }

    rhs[0] += (1/(dx*dx))*2; rhs[y.Size()-1] += (1/(dx*dx))*2;

    d_matrix L(Matrix<double>::Laplacian121(y.Size()));
    L = L * (1/(dx*dx));

    rhs = rhs  - L*y;

    return rhs;
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

    int Psize = 199;
    d_vector y0P(exact_P62_sin(0, Psize));

    d_matrix L(Matrix<double>::Laplacian121(Psize));
    double dx = (1.0-0.0)/(Psize+1);
    L = L * (-1) * (1/(dx*dx));

    TimeIntegration<d_vector> P_62_sin;
    //P_62_sin.Solve(ExponentialRK<d_matrix>::EEuler(L*(-h)), &RHS_P62_sin, h, y0P, t0, t1, steps);
    P_62_sin.Solve(ExponentialRK<d_matrix>::ERK32ZB(L*(-h)), &RHS_P62_sin, h, y0P, t0, t1, steps);
    //P_62_sin.Solve(ClassicalRK::RK4(), &RHS_P62_sin_wL, h, y0P, t0, t1, steps);

    P_62_sin.save_simulation("vector_output.csv", steps);

    return 0;
}
