#include "Matrix.hpp"
#include "Vector.hpp"
#include <mpich/mpi.h>
#include <cassert>
#include <vector>

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
    assert(m>=0);
    assert(n>=0);
    mNumRows = m;
    mNumCols = n;

    mData.resize(mNumRows*mNumCols);
}

int Matrix::NumRows() const
{
    return mNumRows;
}

int Matrix::NumCols() const
{
    return mNumCols;
}

double& Matrix::front()
{
    return mData.front();
}

double& Matrix::operator()(int i, int j)
{
    assert(i > -1);
    assert(i < mNumRows);
    assert(j > -1);
    assert(j < mNumCols);
    return mData[i*mNumCols+j];
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

Vector Matrix::SerialMV(Vector &V)
{
    int ProcNum;        // Number of available processes
    int ProcRank;       // Rank of current process

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    assert(mNumCols == V.Size());
    Vector W(mNumCols);
    double temp;

    if(ProcRank == 0)
    {
        for(int i=0; i<mNumRows; i++)
        {
            temp = 0;
            for(int j=0; j<mNumCols; j++)
            {
                temp += mData[i*mNumCols+j]*V[j];
            }
            W[i]=temp;
        }
    }
    return W;
}

Vector Matrix::operator*(Vector &v) const
{
    int ProcNum;        // Number of available processes
    int ProcRank;       // Rank of current process

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    
    // Number of rows in matrix stripe
    int RowNum = mNumRows/ProcNum;
    // Stripe of the matrix in current process
    Matrix pProcRows(RowNum,mNumCols);
    // Block of result vector in current process
    Vector pProcRes(RowNum);
    // Final resulting vector
    Vector w(mNumRows);

    // Bcast the full vector v1
    MPI_Bcast(&v.front(), v.Size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // Scater the blocks of the matrix
    MPI_Scatter(&mData.front(), RowNum*mNumCols, MPI_DOUBLE,
            &pProcRows.front(), RowNum*mNumCols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Testing if the data distribution of the matrix is correct
    /*
    if(ProcRank == 0)
    {
        std::cout << pProcRows << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(ProcRank == 1)
    {
        std::cout << pProcRows << std::endl;
    }
    */
    for(int i=0; i<ProcNum; i++)
    {
        for(int j=0; j<mNumCols; j++)
        {
            pProcRes[i] += pProcRows(i,j)*v[j];
        }
    }
    // Testing if the parallel multiplication is correct
    /*
    if(ProcRank == 0)
    {
        std::cout << pProcRes << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(ProcRank == 1)
    {
        std::cout << pProcRes << std::endl;
    }
    */
    MPI_Allgather(&pProcRes.front(), RowNum, MPI_DOUBLE,
            &w.front(), RowNum, MPI_DOUBLE, MPI_COMM_WORLD);

    return w;
}
