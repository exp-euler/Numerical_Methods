#ifndef MATRIXHEADERDEF
#define MATRIXHEADERDEF

#include "Vector.hpp"
#include <vector>

class Matrix
{
    private:
        int mNumRows;
        int mNumCols;
        std::vector<std::vector<double>> mData;
    public:
        Matrix(const Matrix& otherMatrix);
        Matrix(int m, int n);
        //~Matrix();
        int NumRows() const;
        int NumCols() const;
        //double& operator()(int i, int j);
        //Matrix& operator=(const Matrix& otherMatrix);
        // assignment
        //Matrix operator+() const;
        //Matrix operator-() const;
        //Matrix operator+(const Matrix& M1) const;
        //Matrix operator-(const Matrix& M1) const;
        double& operator()(int i, int j);
        friend std::ostream& operator<<(std::ostream& output, Matrix& M);
        // multiplication
        //Matrix operator*(double a) const;
        //Matrix operator*(Matrix& M1) const;
        Vector operator*(Vector& v1) const;
};

#endif
