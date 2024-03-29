#ifndef VECTORHEADERDEF
#define VECTORHEADERDEF

#include <cassert>
#include<vector>
#include<iostream>

template<typename DATA_TYPE>
class Vector
{
    private:
        std::vector<DATA_TYPE> mData;
        int mSize;
    public:
        Vector(const Vector<DATA_TYPE>& otherVector);
        Vector(int rows, int cols=1);
        //~Vector();
        int size() const;
        int rows() const;
        int cols() const;
        Vector<DATA_TYPE> cwiseProduct(Vector<DATA_TYPE>& otherVector);
        DATA_TYPE& front();
        DATA_TYPE& operator()(int i);
        Vector<DATA_TYPE>& operator=(const Vector<DATA_TYPE>& otherVector);
        // TODO: Instead of double, use a template here.
        Vector<DATA_TYPE>& operator=(double a);
        Vector<DATA_TYPE> operator+(const Vector<DATA_TYPE>& otherVector);
        Vector<DATA_TYPE> operator*(const Vector<DATA_TYPE>& otherVector);
        Vector<DATA_TYPE> operator-(const Vector<DATA_TYPE>& otherVector);
        Vector<DATA_TYPE> operator+=(const Vector<DATA_TYPE>& otherVector);
        Vector<DATA_TYPE> operator+(const DATA_TYPE scalar);
        Vector<DATA_TYPE> operator*(const DATA_TYPE scalar);

        template<typename D_TYPE>
        friend Vector<D_TYPE> operator+(const D_TYPE scalar, const Vector<D_TYPE>& otherVector);
        template<typename D_TYPE>
        friend Vector<D_TYPE> operator*(const D_TYPE scalar, const Vector<D_TYPE>& otherVector);
        template<typename D_TYPE>
        friend std::ostream& operator<<(std::ostream& output, Vector<D_TYPE>& v);
};


// ###################################################################
// ###################################################################
//                            IMPLEMENTATIONS
// ###################################################################
// ###################################################################

// Overridden copy constructor
template<typename DATA_TYPE>
Vector<DATA_TYPE>::Vector(const Vector<DATA_TYPE>& otherVector)
{
    mSize = otherVector.size();
    // Copying vectors using the std::vector assignment
    mData = otherVector.mData;
}

// Constructor for Vector of a given size
template<typename DATA_TYPE>
Vector<DATA_TYPE>::Vector(int rows, int cols)
{
    assert(rows > 0);
    mSize = rows;

    mData.resize(rows);
}

template<typename DATA_TYPE>
int Vector<DATA_TYPE>::size() const
{
    return mSize;
}

template<typename DATA_TYPE>
int Vector<DATA_TYPE>::rows() const
{
    return mSize;
}

template<typename DATA_TYPE>
int Vector<DATA_TYPE>::cols() const
{
    return 1;
}

template<typename DATA_TYPE>
Vector<DATA_TYPE> Vector<DATA_TYPE>::cwiseProduct(Vector<DATA_TYPE> &otherVector)
{
    assert(mSize == otherVector.mSize);
    Vector<DATA_TYPE> vec(mSize);
    for(int i=0; i<mSize; i++)
    {
        vec.mData[i] = mData[i] * otherVector.mData[i];
    }
    return vec;
}

template<typename DATA_TYPE>
DATA_TYPE& Vector<DATA_TYPE>::front()
{
    return mData.front();
}

template<typename DATA_TYPE>
DATA_TYPE& Vector<DATA_TYPE>::operator()(int i)
{
    assert(i > -1);
    assert(i < mSize);
    return mData[i];
}

template<typename DATA_TYPE>
Vector<DATA_TYPE>& Vector<DATA_TYPE>::operator=(const Vector<DATA_TYPE> &otherVector)
{
    assert(mSize == otherVector.mSize);
    for(int i=0; i<mSize; i++)
    {
        mData[i] = otherVector.mData[i];
    }
    return *this;
}

template<typename DATA_TYPE>
Vector<DATA_TYPE>& Vector<DATA_TYPE>::operator=(double a)
{
    assert(mSize > 0); // TODO: Find better assert?
    for(int i=0; i<mSize; i++)
    {
        mData[i] = a;
    }
    return *this;
}

template<typename DATA_TYPE>
Vector<DATA_TYPE> Vector<DATA_TYPE>::operator+(const Vector<DATA_TYPE> &otherVector)
{
    assert(mSize == otherVector.mSize);
    Vector<DATA_TYPE> vec(mSize);
    for(int i=0; i<mSize; i++)
    {
        vec.mData[i] = mData[i] + otherVector.mData[i];
    }
    return vec;
}

template<typename DATA_TYPE>
Vector<DATA_TYPE> Vector<DATA_TYPE>::operator*(const Vector<DATA_TYPE> &otherVector)
{
    assert(mSize == otherVector.mSize);
    Vector<DATA_TYPE> vec(mSize);
    for(int i=0; i<mSize; i++)
    {
        vec.mData[i] = mData[i] * otherVector.mData[i];
    }
    return vec;
}

template<typename DATA_TYPE>
Vector<DATA_TYPE> Vector<DATA_TYPE>::operator-(const Vector<DATA_TYPE> &otherVector)
{
    assert(mSize == otherVector.mSize);
    Vector<DATA_TYPE> vec(mSize);
    for(int i=0; i<mSize; i++)
    {
        vec.mData[i] = mData[i] - otherVector.mData[i];
    }
    return vec;
}

template<typename DATA_TYPE>
Vector<DATA_TYPE> Vector<DATA_TYPE>::operator+=(const Vector<DATA_TYPE> &otherVector)
{
    assert(mSize == otherVector.mSize);
    for(int i=0; i<mSize; i++)
    {
        mData[i] = mData[i] + otherVector.mData[i];
    }
    return *this;
}

template<typename DATA_TYPE>
Vector<DATA_TYPE> Vector<DATA_TYPE>::operator+(const DATA_TYPE scalar)
{
    Vector<DATA_TYPE> vec(mSize);
    for(int i=0; i<mSize; i++)
    {
        vec.mData[i] = mData[i] + scalar;
    }
    return vec;
}

template<typename DATA_TYPE>
Vector<DATA_TYPE> Vector<DATA_TYPE>::operator*(const DATA_TYPE scalar)
{
    Vector<DATA_TYPE> vec(mSize);
    for(int i=0; i<mSize; i++)
    {
        vec.mData[i] = mData[i] * scalar;
    }
    return vec;
}

template<typename D_TYPE>
Vector<D_TYPE> operator+(const D_TYPE scalar, const Vector<D_TYPE>& otherVector)
{
    Vector<D_TYPE> vec(otherVector.mSize);
    for(int i=0; i<otherVector.mSize; i++)
    {
        vec.mData[i] = otherVector.mData[i] + scalar;
    }
    return vec;
}

template<typename D_TYPE>
Vector<D_TYPE> operator*(const D_TYPE scalar, const Vector<D_TYPE>& otherVector)
{
    Vector<D_TYPE> vec(otherVector.mSize);
    for(int i=0; i<otherVector.mSize; i++)
    {
        vec.mData[i] = otherVector.mData[i] * scalar;
    }
    return vec;
}

template<typename D_TYPE>
std::ostream& operator<<(std::ostream& output, Vector<D_TYPE>& v)
{
    // Format for when outputting to console
    if(&output == &std::cout){
        output << std::endl;
        for(int i=0; i<v.size(); i++)
        {
            std::cout << v(i) << std::endl;
        }
        return output; // return std::ostream so that we can have chain call of <<
    }

    // Format for when outputting to csv file
    for(int i=0; i<v.size(); i++)
    {
        output << v(i) << ",";
    }
    return output; // return std::ostream so that we can have chain call of <<
}

#endif
