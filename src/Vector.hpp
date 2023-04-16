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
        Vector(int size);
        //~Vector();
        int Size() const;
        DATA_TYPE& front();
        DATA_TYPE& operator[](int i);
        Vector<DATA_TYPE>& operator=(const Vector<DATA_TYPE>& otherVector);
        Vector<DATA_TYPE>& operator+(const Vector<DATA_TYPE>& otherVector);
        Vector<DATA_TYPE>& operator+=(const Vector<DATA_TYPE>& otherVector);
        Vector<DATA_TYPE>& operator+(const DATA_TYPE scalar);
        Vector<DATA_TYPE>& operator*(const DATA_TYPE scalar);

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
    mSize = otherVector.Size();
    // Copying vectors using the std::vector assignment
    mData = otherVector.mData;
}

// Constructor for Vector of a given size
template<typename DATA_TYPE>
Vector<DATA_TYPE>::Vector(int size)
{
    assert(size > 0);
    mSize = size;

    mData.resize(size);
}

template<typename DATA_TYPE>
int Vector<DATA_TYPE>::Size() const
{
    return mSize;
}

template<typename DATA_TYPE>
DATA_TYPE& Vector<DATA_TYPE>::front()
{
    return mData.front();
}

template<typename DATA_TYPE>
DATA_TYPE& Vector<DATA_TYPE>::operator[](int i)
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
Vector<DATA_TYPE>& Vector<DATA_TYPE>::operator+(const Vector<DATA_TYPE> &otherVector)
{
    assert(mSize == otherVector.mSize);
    for(int i=0; i<mSize; i++)
    {
        mData[i] = mData[i] + otherVector.mData[i];
    }
    return *this;
}

template<typename DATA_TYPE>
Vector<DATA_TYPE>& Vector<DATA_TYPE>::operator+=(const Vector<DATA_TYPE> &otherVector)
{
    assert(mSize == otherVector.mSize);
    for(int i=0; i<mSize; i++)
    {
        mData[i] = mData[i] + otherVector.mData[i];
    }
    return *this;
}

template<typename DATA_TYPE>
Vector<DATA_TYPE>& Vector<DATA_TYPE>::operator+(const DATA_TYPE scalar)
{
    for(int i=0; i<mSize; i++)
    {
        mData[i] = mData[i] + scalar;
    }
    return *this;
}

template<typename DATA_TYPE>
Vector<DATA_TYPE>& Vector<DATA_TYPE>::operator*(const DATA_TYPE scalar)
{
    for(int i=0; i<mSize; i++)
    {
        mData[i] = mData[i] * scalar;
    }
    return *this;
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
    output << std::endl;
    for(int i=0; i<v.Size(); i++)
    {
        std::cout << v[i] << std::endl;
    }
    return output; // return std::ostream so that we can have chain call of <<
}

#endif
