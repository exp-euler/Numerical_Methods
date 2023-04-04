#include "Matrix.hpp"
#include "Vector.hpp"

int main(){
    Matrix<double> A(3,3);
    Matrix<double> B(3,3);
    Matrix<double> C(3,3);

    A(0,0)=1; A(0,1)=2; A(0,2)=3;
    A(1,0)=1; A(1,1)=2; A(1,2)=3;
    A(2,0)=1; A(2,1)=2; A(2,2)=3;

    B(0,0)=1; B(0,1)=0; B(0,2)=0;
    B(1,0)=0; B(1,1)=1; B(1,2)=0;
    B(2,0)=0; B(2,1)=0; B(2,2)=1;

    C = A*B;

    std::cout << C << "\n";

    return 0;
}
