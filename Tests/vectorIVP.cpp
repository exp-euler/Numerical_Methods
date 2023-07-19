//TODO: Add the exact solution and do a check
#include "TimeIntegration.hpp"
#include "Vector.hpp"

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

int main() {
    TimeIntegration<d_vector> system;
    int steps = 1000;
    double t0 = 0;
    double t1 = 1;
    double h = (t1-t0)/steps;
    d_vector y0(2);
    y0[0] = 0; y0[1] = -4;
    system.Solve(ClassicalRK::Euler(), &RHS, h, y0, t0, t1, steps);

    d_vector error(exact(t1)-system.getY()[system.getY().size()-1]);

    std::cout << "Error is: " << error << "\n";
    return 0;
}
