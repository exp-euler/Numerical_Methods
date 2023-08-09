#include "MatrixFunctions.hpp"

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


int main(){
    d_matrix A(3,3);
    d_matrix B(3,3);
    d_matrix C(3,3);
    d_matrix L(2,2);

    A(0,0)=1; A(0,1)=2; A(0,2)=3;
    A(1,0)=1; A(1,1)=2; A(1,2)=3;
    A(2,0)=1; A(2,1)=2; A(2,2)=-30;

    B(0,0)=1; B(0,1)=0; B(0,2)=0;
    B(1,0)=0; B(1,1)=1; B(1,2)=0;
    B(2,0)=7; B(2,1)=0; B(2,2)=1;

    L(0,0)=-2; L(0,1)=1;
    L(1,0)=998; L(1,1)=-999;

    C = A*B;
    C = A+B;

    d_matrix phi0(2,2);
    d_matrix phi1(2,2);
    d_matrix phi2(2,2);
    d_matrix phi3(2,2);
    std::vector<d_matrix> phi_f = {phi0, phi1, phi2, phi3};
    LinearAlgebra::phi_functions(phi_f, 3, L*0.001);

    //Matrix<double> Lap(Matrix<double>::Laplacian121(10));
    //std::cout << Lap << "\n";
    std::cout << "##########################" << "\n";
    std::cout << phi_f[0] << "\n";
    std::cout << phi_f[1] << "\n";
    std::cout << phi_f[2] << "\n";
    std::cout << phi_f[3] << "\n";

    return 0;
}
