#include "TimeIntegration.hpp"
#include "Vector.hpp"
#include <cmath>
#include <getopt.h>

typedef Vector<double> d_vector;

d_vector RHS(double t, d_vector y) {
    d_vector rhs(y.Size());
    rhs[0] = y[0] + 2*y[1];
    rhs[1] = 3*y[0] + 2*y[1];
    return rhs;
}

d_vector exact(double t) {
    d_vector sol(2);
    sol[0] = (-8.0/5.0)*std::exp(-t)*(-1) - (4.0/5.0)*std::exp(4*t)*2;
    sol[1] = (-8.0/5.0)*std::exp(-t)*1 - (4.0/5.0)*std::exp(4*t)*3;
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
    y0[0] = 0; y0[1] = -4;

    TimeIntegration<d_vector> system;
    system.Solve(ClassicalRK::Euler(), &RHS, h, y0, t0, t1, steps);

    d_vector error(exact(t1)-system.getY()[system.getY().size()-1]);
    std::cout << "Error is: " << error << "\n";

    system.save_simulation("vector_output.csv", steps);
    return 0;
}
