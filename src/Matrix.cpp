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
    assert(m>0);
    assert(n>0);
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

Vector Matrix::SerialMV(Vector &V)
{
    assert(mNumCols == V.Size());
    Vector W(mNumRows);
    double temp;

    for(int i=0; i<mNumRows; i++)
    {
        temp = 0;
        for(int j=0; j<mNumCols; j++)
        {
            temp += mData[i*mNumCols+j]*V[j];
        }
        W[i]=temp;
    }
    return W;
}

Matrix Matrix::SerialMM(Matrix &M)
{
    assert(mNumCols == M.mNumRows);
    Matrix N(mNumRows, M.mNumCols);
    double temp;

    for(int i=0; i<N.mNumRows; i++)
    {
        for(int j=0; j<N.mNumCols; j++)
        {
            N(i,j)=0;
            for(int k=0; k<mNumCols; k++)
            {
                N(i,j) += mData[i*mNumCols+k]*M(k,j);
            }
        }
    }
    return N;
}

double& Matrix::operator()(int i, int j)
{
    assert(i > -1);
    assert(i < mNumRows);
    assert(j > -1);
    assert(j < mNumCols);
    return mData[i*mNumCols+j];
}

Matrix& Matrix::operator=(const Matrix &otherMatrix)
{
    assert(mNumCols == otherMatrix.mNumCols);
    assert(mNumRows == otherMatrix.mNumRows);
    for(int i=0; i<mNumRows*mNumCols; i++)
    {
        mData[i] = otherMatrix.mData[i];
    }
    return *this;
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

Vector Matrix::operator*(Vector &v) const
{
    int ProcNum;        // Number of available processes
    int ProcRank;       // Rank of current process

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    
    // Number of rows that haven't been distributed yet
    int RestRows = mNumRows;
    for(int i=0; i<ProcRank; i++)
    {
        RestRows = RestRows - RestRows/(ProcNum-i);
    }
    // Number of rows in matrix stripe for best load balancing
    int RowNum = RestRows/(ProcNum-ProcRank);
    // Stripe of the matrix in current process
    Matrix ProcRows(RowNum,mNumCols);
    // Block of result vector in current process
    Vector ProcRes(RowNum);
    // Final resulting vector
    Vector w(mNumRows);

    // Bcast the full vector v1
    MPI_Bcast(&v.front(), v.Size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*         Load balancing of scattering of matrix M             */
    int *pSendInd; // Index of first data element sent to a proc via Scatterv
    int *pSendNum; // Num of elements sent to a proc via Scatterv
    RestRows = mNumRows;
    // Alloc memory for temp objects for scattering
    pSendInd = new int [ProcNum];
    pSendNum = new int [ProcNum];
    // Determine how the matrix is distributed to processes.
    RowNum = mNumRows/ProcNum;
    pSendInd[0] = 0;
    pSendNum[0] = RowNum*mNumCols;
    for(int i=1; i<ProcNum; i++)
    {
        RestRows -= RowNum;
        RowNum = RestRows/(ProcNum-i);
        pSendNum[i] = RowNum*mNumCols;
        pSendInd[i] = pSendInd[i-1]+pSendNum[i-1];
    }

    // Scater the blocks of the matrix
    MPI_Scatterv(&mData.front(), pSendNum, pSendInd, MPI_DOUBLE,
            &ProcRows.front(), pSendNum[ProcRank], MPI_DOUBLE,
            0, MPI_COMM_WORLD);
    delete [] pSendNum;
    delete [] pSendInd;

    // Testing if the data distribution of the matrix is correct
    /*
    if(ProcRank == 0)
    {
        std::cout << ProcRows << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(ProcRank == 1)
    {
        std::cout << ProcRows << std::endl;
    }
    */
    for(int i=0; i<ProcRes.Size(); i++)
    {
        for(int j=0; j<mNumCols; j++)
        {
            ProcRes[i] += ProcRows(i,j)*v[j];
        }
    }
    // Testing if the parallel multiplication is correct
    /*
    if(ProcRank == 0)
    {
        std::cout << ProcRes << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(ProcRank == 1)
    {
        std::cout << ProcRes << std::endl;
    }
    */

    /*         Load balancing of gathering the resulting vector w             */
    int *pReceiveInd; // Index of first data element from current process to w
    int *pReceiveNum; // Numb of elements current process sends with Allgatherv
    // Alloc memory for temp objects for scattering
    pReceiveInd = new int [ProcNum];
    pReceiveNum = new int [ProcNum];
    // Determine how the matrix is distributed to processes.
    RestRows = mNumRows;
    pReceiveInd[0] = 0;
    pReceiveNum[0] = mNumRows/ProcNum;
    for(int i=1; i<ProcNum; i++)
    {
        RestRows -= pReceiveNum[i-1];
        pReceiveNum[i] = RestRows/(ProcNum - i);
        pReceiveInd[i] = pReceiveInd[i-1] + pReceiveNum[i-1];
    }

    MPI_Allgatherv(&ProcRes.front(), pReceiveNum[ProcRank], MPI_DOUBLE,
            &w.front(), pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
    
    delete [] pReceiveNum;
    delete [] pReceiveInd;

    return w;
}

Matrix Matrix::operator*(Matrix &M) const
{
    int ProcNum;        // Number of available processes
    int ProcRank;       // Rank of current process

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    // Number of rows of first matrix that haven't been distributed yet
    int RestRows = mNumRows;
    for(int i=0; i<ProcRank; i++)
    {
        RestRows = RestRows - RestRows/(ProcNum-i);
    }
    // Number of rows in matrix stripe for best load balancing
    int RowNum = RestRows/(ProcNum-ProcRank);
    // Stripe of the matrix in current process
    Matrix ProcRows(RowNum,mNumCols);
    // Block of result matrix in current process
    Matrix ProcRes(RowNum,M.mNumCols);
    // Final resulting matrix
    Matrix N(mNumRows,M.mNumCols);


    // Bcast the full matrix M
    MPI_Bcast(&M.front(), M.mNumRows*M.mNumCols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*         Load balancing of scattering of matrix M             */
    int *pSendInd; // Index of first data element sent to a proc via Scatterv
    int *pSendNum; // Num of elements sent to a proc via Scatterv
    RestRows = mNumRows;
    // Alloc memory for temp objects for scattering
    pSendInd = new int [ProcNum];
    pSendNum = new int [ProcNum];
    // Determine how the matrix is distributed to processes.
    RowNum = mNumRows/ProcNum;
    pSendInd[0] = 0;
    pSendNum[0] = RowNum*mNumCols;
    for(int i=1; i<ProcNum; i++)
    {
        RestRows -= RowNum;
        RowNum = RestRows/(ProcNum-i);
        pSendNum[i] = RowNum*mNumCols;
        pSendInd[i] = pSendInd[i-1]+pSendNum[i-1];
    }

    // Scater the blocks of the matrix
    MPI_Scatterv(&mData.front(), pSendNum, pSendInd, MPI_DOUBLE,
            &ProcRows.front(), pSendNum[ProcRank], MPI_DOUBLE,
            0, MPI_COMM_WORLD);
    delete [] pSendNum;
    delete [] pSendInd;

    // Testing if the data distribution of the matrix is correct
    /*
    if(ProcRank == 0)
    {
        std::cout << ProcRows << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(ProcRank == 1)
    {
        std::cout << ProcRows << std::endl;
    }
    */
    for(int i=0; i<ProcRes.mNumRows; i++)
    {
        for(int j=0; j<ProcRes.mNumCols; j++)
        {
            ProcRes(i,j)=0;
            for(int k=0; k<mNumCols; k++)
            {
                ProcRes(i,j) += ProcRows(i,k)*M(k,j);
            }
        }
    }
    // Testing if the parallel multiplication is correct
    /*
    if(ProcRank == 0)
    {
        std::cout << ProcRes << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(ProcRank == 1)
    {
        std::cout << ProcRes << std::endl;
    }
    */

    /*         Load balancing of gathering the resulting matrix N             */
    int *pReceiveInd; // Index of first data element from current process to N
    int *pReceiveNum; // Numb of elements current process sends with Allgatherv
    // Alloc memory for temp objects for scattering
    pReceiveInd = new int [ProcNum];
    pReceiveNum = new int [ProcNum];
    // Determine how the matrix is distributed to processes.
    RestRows = mNumRows;
    pReceiveInd[0] = 0;
    pReceiveNum[0] = (RestRows/ProcNum)*M.mNumCols;
    for(int i=1; i<ProcNum; i++)
    {
        RestRows -= pReceiveNum[i-1]/M.mNumCols;
        pReceiveNum[i] = (RestRows/(ProcNum - i))*M.mNumCols;
        pReceiveInd[i] = pReceiveInd[i-1] + pReceiveNum[i-1];
    }

    MPI_Allgatherv(&ProcRes.front(), pReceiveNum[ProcRank], MPI_DOUBLE,
            &N.front(), pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
    
    delete [] pReceiveNum;
    delete [] pReceiveInd;

    return N;
}
