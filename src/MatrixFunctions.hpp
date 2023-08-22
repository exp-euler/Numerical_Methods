#ifndef MATRIXFUNCTIONS
#define MATRIXFUNCTIONS
#include <vector>
#include <cassert>
#include <bits/stdc++.h>
#ifdef EIGEN_YES

#include <Eigen/Dense>
typedef Eigen::MatrixXd matrix;
typedef Eigen::VectorXd vector;

#else

#include "Matrix.hpp"
#include "Vector.hpp"
typedef Matrix<double> matrix;
typedef Vector<double> vector;

#endif


namespace LinearAlgebra {
    void Taylor(matrix &T, matrix &L, int k);

    int scale(matrix &A);

    void unscale(std::vector<matrix> &phi, matrix L, int j, int l);

    void phi_functions(std::vector<matrix> &phi, int l, matrix A);

    void phi_functions(std::vector<vector> &phi, int l, vector A);
}
#endif //MATRIXFUNCTIONS
