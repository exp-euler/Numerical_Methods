#ifndef MATRIXFUNCTIONS
#define MATRIXFUNCTIONS
#include "Matrix.hpp"
#include <vector>
#include <cassert>
#include <bits/stdc++.h>

typedef Matrix<double> matrix;
namespace LinearAlgebra {
    void Taylor(matrix &T, matrix &L, int k);

    int scale(matrix &A);

    void unscale(std::vector<matrix> &phi, matrix L, int j, int l);

    void phi_functions(std::vector<matrix> &phi, int l, matrix A);

}
#endif //MATRIXFUNCTIONS
