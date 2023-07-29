#ifndef TABLEAUS
#define TABLEAUS

#include "Matrix.hpp"
#include "MatrixFunctions.hpp"
#include<vector>

// Provides a way to store an RK method by its Butcher Tableau
// Load the c_i (steps), a_ij (weights) and b_i (final weights).
class ClassicalRK
{
    // All class methods are public for ease of working.
    // Always pass a Tableau Object by const reference to avoid editing tableaus.
    public:
    // Add more RK methods as named constructors like below
    static ClassicalRK Euler();
    static ClassicalRK RK4();

    std::vector<double> c;
    // e is a vector with transformation factors. Len = Len(c) + 1
    std::vector<double> e; // In case of RK it is 1
    std::vector<std::vector<double>> a;
    std::vector<double> b;

    ClassicalRK(std::vector<double> steps,
                std::vector<double> t_factors,
                std::vector<std::vector<double>> weights,
                std::vector<double> weightsHigher);
};

template<typename DATA_TYPE>
class ExponentialRK
{
    public:
    // Add more ERK methods as named constructors like below
    static ExponentialRK EEuler(DATA_TYPE L);
    static ExponentialRK ERK32ZB(DATA_TYPE L);

    std::vector<double> c;
    // e is a vector with transformation factors.
    std::vector<DATA_TYPE> e; // In case of ERK it is exp(-c*L*h)
    std::vector<std::vector<DATA_TYPE>> a;
    std::vector<DATA_TYPE> b;

    ExponentialRK(DATA_TYPE L,
                  std::vector<double> steps,
                  std::vector<DATA_TYPE> t_factors,
                  std::vector<std::vector<DATA_TYPE>> weights,
                  std::vector<DATA_TYPE> weightsHigher);
};

template<typename DATA_TYPE>
ExponentialRK<DATA_TYPE>::ExponentialRK( DATA_TYPE L,
             std::vector<double> steps,
             std::vector<DATA_TYPE> t_factors,
             std::vector<std::vector<DATA_TYPE>> weights,
             std::vector<DATA_TYPE> weightsHigher)
{
    c.reserve(steps.size());
    e.reserve(t_factors.size());
    a.reserve(weights.size());
    b.reserve(weightsHigher.size());

    c = steps;
    e = t_factors;
    a = weights;
    b = weightsHigher;
}

template<typename DATA_TYPE>
ExponentialRK<DATA_TYPE> ExponentialRK<DATA_TYPE>::EEuler(DATA_TYPE L)
{
    int size_L = L.NumRows();
    Matrix<double> phi0(size_L,size_L);
    Matrix<double> phi1(size_L,size_L);
    std::vector<DATA_TYPE> phi = {phi0, phi1};

    LinearAlgebra::phi_functions(phi, 1, L);

    return ExponentialRK(L, {},{phi[0]},{},{phi[1]});
}

template<typename DATA_TYPE>
ExponentialRK<DATA_TYPE> ExponentialRK<DATA_TYPE>::ERK32ZB(DATA_TYPE L)
{
    int size_L = L.NumRows();
    Matrix<double> phi00(size_L,size_L);
    Matrix<double> phi10(size_L,size_L);
    Matrix<double> phi20(size_L,size_L);
    Matrix<double> phi01(size_L,size_L);
    Matrix<double> phi11(size_L,size_L);
    Matrix<double> phi21(size_L,size_L);
    Matrix<double> phi0(size_L,size_L);
    Matrix<double> phi1(size_L,size_L);
    Matrix<double> phi2(size_L,size_L);
    Matrix<double> phi3(size_L,size_L);

    std::vector<DATA_TYPE> phi = {phi0, phi1, phi2, phi3};
    std::vector<DATA_TYPE> phi_c0 = {phi00, phi10, phi20};
    std::vector<DATA_TYPE> phi_c1 = {phi01, phi11, phi21};

    LinearAlgebra::phi_functions(phi, 3, L);
    LinearAlgebra::phi_functions(phi_c0, 2, L*0.5);
    LinearAlgebra::phi_functions(phi_c1, 2, L*0.75);

    return ExponentialRK(L,
        {0.5, 0.75},
        {phi_c0[0], phi_c1[0], phi[0]},
        {{phi_c0[1]*0.5}, {phi_c1[1]*0.75 - phi_c1[2]*(9.0/8.0) - phi_c0[2]*(3.0/8.0), phi_c1[2]*(9.0/8.0) + phi_c0[2]*(3.0/8.0)}},
        {phi[1] - phi[2]*0.75 + phi[3]*0.25 - phi[2]*(5.0/6.0) - phi[3]*(1.0/6), phi[2]*0.75 - phi[3]*0.25, phi[2]*(5.0/6.0) + phi[3]*(1.0/6)});
}



#endif //TABLEAUS
