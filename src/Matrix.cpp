#include "Matrix.hpp"

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
