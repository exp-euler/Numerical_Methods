#include "Matrix.hpp"
#include "Vector.hpp"
#include <cassert>

// Overridden copy constructor
Matrix::Matrix(const Matrix& otherMatrix)
{
    mNumRows = otherMatrix.mNumRows;
    mNumCols = otherMatrix.mNumCols;

    // Copying matrixes using std::vector assignment
    mData = otherMatrix.mData;
}

// Constructor for Matrix of a given size
Matrix::Matrix(int m, int n)
{
    assert(m>0);
    assert(n>0);
    mNumRows = m;
    mNumCols = n;

    mData.reserve(mNumRows);

    for(int i =0; i<mNumRows; i++)
    {
        mData[i].reserve(mNumCols);
        for(int j=0; j<mNumCols; j++)
        {
            mData[i][j] = 0.0;
        }
    }
}

int Matrix::NumRows() const
{
    return mNumRows;
}

int Matrix::NumCols() const
{
    return mNumCols;
}

double& Matrix::operator()(int i, int j)
{
    assert(i > -1);
    assert(i < mNumRows);
    assert(j > -1);
    assert(j < mNumCols);
    return mData[i][j];
}

std::ostream& operator<<(std::ostream& output, Matrix& M)
{
    output << std::endl;
    for(int i=0; i<M.NumRows(); i++)
    {
        for(int j=0; j<M.NumCols();j++)
        {
            std::cout << M(i,j) << " ";
        }
        std::cout << std::endl;
    }
    return output; // return std::ostream so that we can have chain call of <<
}

Vector Matrix::operator*(Vector &v1) const
{
    assert(mNumCols == v1.Size());
    Vector V(mNumCols);
    double temp;
    for(int i=0; i<mNumRows; i++)
    {
        temp = 0;
        for(int j=0; j<mNumCols; j++)
        {
            temp += mData[i][j]*v1[j];
        }
        V[i]=temp;
    }
    return V;
}
