//TODO: Add the exact solution and do a check
#include "TimeIntegration.hpp"
#include "Vector.hpp"

typedef Vector<double> d_vector;

d_vector RHS(double t, d_vector y) {
    return y + 1;
}

int main() {
    TimeIntegration<d_vector> system;
    int steps = 50;
    double t0 = 0;
    double t1 = 1;
    double h = (t1-t0)/steps;
    d_vector y0(2);
    y0[0] = 1; y0[1] = 1;
    system.Solve(ClassicalRK::Euler(), &RHS, h, y0, t0, t1, steps);

    // Test
    d_vector y_test(2);
    y_test = y0;
    y_test += y0;
    std::cout << y_test << "\n";
    std::cout << y0 << "\n";
    return 0;
}
