#include "Matrix.hpp"
#include "Vector.hpp"
#include "MatrixFunctions.hpp"

int main(){
    Matrix<double> A(3,3);
    Matrix<double> B(3,3);
    Matrix<double> C(3,3);

    A(0,0)=1; A(0,1)=2; A(0,2)=3;
    A(1,0)=1; A(1,1)=2; A(1,2)=3;
    A(2,0)=1; A(2,1)=2; A(2,2)=-30;

    B(0,0)=1; B(0,1)=0; B(0,2)=0;
    B(1,0)=0; B(1,1)=1; B(1,2)=0;
    B(2,0)=7; B(2,1)=0; B(2,2)=1;

    C = A*B;
    C = A+B;

    Matrix<double> phi1(3,3);
    Matrix<double> phi2(3,3);
    std::vector<Matrix<double>> phi_f = {phi1, phi2};
    LinearAlgebra::phi_functions(phi_f, 1, A);

    std::cout << A << "\n";
    std::cout << "##########################" << "\n";
    std::cout << phi_f[1] << "\n";

    return 0;
}
