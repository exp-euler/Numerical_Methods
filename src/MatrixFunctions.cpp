#include "MatrixFunctions.hpp"

typedef Matrix<double> matrix;
namespace LinearAlgebra {
    void Taylor(matrix &T, matrix &L, int k) {
        std::vector<double> coeff = {1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0, 1.0/120.0, 1.0/720.0, 1.0/5040.0, 1.0/40320.0, 1.0/362880.0, 1.0/3628800.0};
        matrix temp(matrix::Diagonal(L.NumRows(), L.NumCols(), 1.0));
        T = temp*coeff[0];
        for(int i=1; i<=k; i++) {
            temp = temp*L;
            T = T + temp*coeff[i];
        }
    }

    int scale(matrix &A) {
    double normA = A.inf_norm();
    int j = normA > 0 ? 1+round(log(normA/4)/log(2)) : 0;
    return j >= 0 ? j : 0;
    }

    void unscale(std::vector<matrix> &phi, matrix L, int j, int l) {
    for(int m=0; m <= l; ++m)
        Taylor(phi[m],L,10);
    
    // Unscaling relations taken form:
    // http://www.math.ntnu.no/preprint/numerics/2005/N4-2005.pdf
    for(int k=1; k <= j; ++k) {
        if(l >= 3)
            phi[3] = (phi[1]*phi[2]+phi[3]*2+phi[2])*0.125;
        if(l >= 2)
            phi[2] = (phi[1]*phi[1]+phi[2]*2)*0.25;
        if(l >= 1)
            phi[1] = (phi[0]*phi[1]+phi[1])*0.5;
        phi[0] = phi[0]*phi[0];
        }

    }

    void phi_functions(std::vector<matrix> &phi, int l, matrix A) {
        int j = scale(A);
        j = 40; // Problem in accuracy is the scaling... Must be high enough?
        unscale(phi, A*(1.0/double(pow(2,j))), j, l);
    }

}
