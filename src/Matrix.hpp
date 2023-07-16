#ifndef MATRIXHEADERDEF
#define MATRIXHEADERDEF

#include "Vector.hpp"
#include <vector>
#include <bits/stdc++.h>
#include <cassert>

#ifdef HYBRID_CPU
#include <mpi.h>
#include <omp.h>
#endif

template<typename DATA_TYPE>
class Matrix
{
    private:
        int mNumRows;
        int mNumCols;
        std::vector<DATA_TYPE> mData;
    public:
        Matrix(const Matrix& otherMatrix);
        Matrix(int m, int n);
        static Matrix<DATA_TYPE> Diagonal(int m, int n, DATA_TYPE k);
        //~Matrix();
        int NumRows() const;
        int NumCols() const;

        // The following is used in the MPI part of this file
        // TODO: Do this step in a different way
        DATA_TYPE& front();
        Matrix<DATA_TYPE> operator*(Matrix<DATA_TYPE>& M) const;
        Vector<DATA_TYPE> operator*(Vector<DATA_TYPE>& v1) const;

        DATA_TYPE inf_norm() const;

        Vector<DATA_TYPE> SerialMV(Vector<DATA_TYPE>& V);
        Matrix<DATA_TYPE> SerialMM(Matrix<DATA_TYPE>& N);

        // Addition/Substraction
        //Matrix operator-() const;
        Matrix<DATA_TYPE> operator+(const Matrix<DATA_TYPE>& M) const;
        Matrix<DATA_TYPE> operator-(const Matrix<DATA_TYPE>& M) const;

        // multiplication
        Matrix<DATA_TYPE> operator*(DATA_TYPE a) const;
        Matrix<DATA_TYPE> operator*(const Matrix<DATA_TYPE>& M) const;
        Vector<DATA_TYPE> operator*(const Vector<DATA_TYPE>& v1) const;

        // Overload both const and non-const version so we can use it
        // on a const this or equivalent referance.
        DATA_TYPE& operator()(int i, int j);
        const DATA_TYPE& operator()(int i, int j) const;

        Matrix<DATA_TYPE>& operator=(const Matrix<DATA_TYPE> &otherMatrix);
        template<typename D_TYPE>
        friend std::ostream& operator<<(std::ostream& output, Matrix<D_TYPE>& M);
};


// ###################################################################
// ###################################################################
//                            IMPLEMENTATIONS
// ###################################################################
// ###################################################################

// Overridden copy constructor
template<typename DATA_TYPE>
Matrix<DATA_TYPE>::Matrix(const Matrix<DATA_TYPE>& otherMatrix)
{
    mNumRows = otherMatrix.mNumRows;
    mNumCols = otherMatrix.mNumCols;

    // Copying matrixes using std::vector assignment
    mData = otherMatrix.mData;
}

// Constructor for Matrix of a given size
template<typename DATA_TYPE>
Matrix<DATA_TYPE>::Matrix(int m, int n)
{
    assert(m>0);
    assert(n>0);
    mNumRows = m;
    mNumCols = n;

    mData.resize(mNumRows*mNumCols);
}

// Named constructor to give identity matrix
template<typename DATA_TYPE>
Matrix<DATA_TYPE> Matrix<DATA_TYPE>::Diagonal(int m, int n, DATA_TYPE k)
{
    assert(m==n);
    Matrix<DATA_TYPE> D(m,n);
    for(int i=0; i<m; i++)
        D(i,i) = k;
    return D;
}

template<typename DATA_TYPE>
int Matrix<DATA_TYPE>::NumRows() const
{
    return mNumRows;
}

template<typename DATA_TYPE>
int Matrix<DATA_TYPE>::NumCols() const
{
    return mNumCols;
}

template<typename DATA_TYPE>
DATA_TYPE& Matrix<DATA_TYPE>::front()
{
    return mData.front();
}

template<typename DATA_TYPE>
DATA_TYPE Matrix<DATA_TYPE>::inf_norm() const
{
    DATA_TYPE row_sums[mNumRows] = {};
    for(int i=0; i<mNumRows; i++)
    {
        for(int j=0; j<mNumCols; j++)
            row_sums[i] += std::abs(mData[i*mNumCols+j]);
    }

    return *std::max_element(row_sums, row_sums+mNumRows);
}

template<typename DATA_TYPE>
Vector<DATA_TYPE> Matrix<DATA_TYPE>::SerialMV(Vector<DATA_TYPE> &V)
{
    assert(mNumCols == V.Size());
    Vector<DATA_TYPE> W(mNumRows);
    DATA_TYPE temp;

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

template<typename DATA_TYPE>
Matrix<DATA_TYPE> Matrix<DATA_TYPE>::SerialMM(Matrix<DATA_TYPE> &M)
{
    assert(mNumCols == M.mNumRows);
    Matrix<DATA_TYPE> N(mNumRows, M.mNumCols);

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

template<typename DATA_TYPE>
const DATA_TYPE& Matrix<DATA_TYPE>::operator()(int i, int j) const
{
    assert(i > -1);
    assert(i < mNumRows);
    assert(j > -1);
    assert(j < mNumCols);
    return mData[i*mNumCols+j];
}

template<typename DATA_TYPE>
DATA_TYPE& Matrix<DATA_TYPE>::operator()(int i, int j)
{
    assert(i > -1);
    assert(i < mNumRows);
    assert(j > -1);
    assert(j < mNumCols);
    return mData[i*mNumCols+j];
}

template<typename DATA_TYPE>
Matrix<DATA_TYPE>& Matrix<DATA_TYPE>::operator=(const Matrix<DATA_TYPE> &otherMatrix)
{
    assert(mNumCols == otherMatrix.mNumCols);
    assert(mNumRows == otherMatrix.mNumRows);
    for(int i=0; i<mNumRows*mNumCols; i++)
    {
        mData[i] = otherMatrix.mData[i];
    }
    return *this;
}

template<typename D_TYPE>
std::ostream& operator<<(std::ostream& output, Matrix<D_TYPE>& M)
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

#ifdef SERIAL

template<typename DATA_TYPE>
Matrix<DATA_TYPE> Matrix<DATA_TYPE>::operator*(DATA_TYPE a) const
{
    Matrix<DATA_TYPE> N(mNumRows, mNumCols);

    for(int i=0; i<mNumRows; i++)
    {
        for(int j=0; j<mNumCols; j++)
        {
            N(i,j)=a*mData[i*mNumCols+j];
        }
    }
    return N;
}


template<typename DATA_TYPE>
Vector<DATA_TYPE> Matrix<DATA_TYPE>::operator*(const Vector<DATA_TYPE> &V) const
{
    assert(mNumCols == V.Size());
    Vector<DATA_TYPE> W(mNumRows);
    DATA_TYPE temp;

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

template<typename DATA_TYPE>
Matrix<DATA_TYPE> Matrix<DATA_TYPE>::operator*(const Matrix<DATA_TYPE> &M) const
{
    assert(mNumCols == M.mNumRows);
    Matrix<DATA_TYPE> N(mNumRows, M.mNumCols);

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

template<typename DATA_TYPE>
Matrix<DATA_TYPE> Matrix<DATA_TYPE>::operator+(const Matrix<DATA_TYPE> &M) const
{
    assert(mNumCols == M.mNumRows);
    Matrix<DATA_TYPE> N(mNumRows, M.mNumCols);

    for(int i=0; i<N.mNumRows; i++)
    {
        for(int j=0; j<N.mNumCols; j++)
        {
            N(i,j)=mData[i*mNumCols+j] + M(i,j);
        }
    }
    return N;
}

template<typename DATA_TYPE>
Matrix<DATA_TYPE> Matrix<DATA_TYPE>::operator-(const Matrix<DATA_TYPE> &M) const
{
    assert(mNumCols == M.mNumRows);
    Matrix<DATA_TYPE> N(mNumRows, M.mNumCols);

    for(int i=0; i<N.mNumRows; i++)
    {
        for(int j=0; j<N.mNumCols; j++)
        {
            N(i,j)=mData[i*mNumCols+j] - M(i,j);
        }
    }
    return N;
}


#endif

#ifdef HYBRID_CPU
template<typename DATA_TYPE>
Vector<DATA_TYPE> Matrix<DATA_TYPE>::operator*(Vector<DATA_TYPE> &v) const
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
    Matrix<DATA_TYPE> ProcRows(RowNum,mNumCols);
    // Block of result vector in current process
    Vector<DATA_TYPE> ProcRes(RowNum);
    // Final resulting vector
    Vector<DATA_TYPE> w(mNumRows);

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

template<typename DATA_TYPE>
Matrix<DATA_TYPE> Matrix<DATA_TYPE>::operator*(Matrix<DATA_TYPE> &M) const
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
    Matrix<DATA_TYPE> ProcRows(RowNum,mNumCols);
    // Block of result matrix in current process
    Matrix<DATA_TYPE> ProcRes(RowNum,M.mNumCols);
    // Final resulting matrix
    Matrix<DATA_TYPE> N(mNumRows,M.mNumCols);


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
    
    // Each process runs in parallel using multithreading
    #pragma omp parallel num_threads(omp_get_max_threads()/ProcNum)
    for(int i=0; i<ProcRes.mNumRows; i++)
    {
        #pragma omp for
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
#endif

#endif
