#include "Vector.hpp"
#include <cassert>

// Overridden copy constructor
Vector::Vector(const Vector& otherVector)
{
    mSize = otherVector.Size();
    // Copying vectors using the std::vector assignment
    mData = otherVector.mData;
}

// Constructor for Vector of a given size
Vector::Vector(int size)
{
    assert(size > 0);
    mSize = size;
    mData.reserve(size);
    for(int i=0; i<mSize; i++)
    {
        mData[i] = 0.0;
    }
}

int Vector::Size() const
{
    return mSize;
}

double& Vector::operator[](int i)
{
    assert(i > -1);
    assert(i < mSize);
    return mData[i];
}

Vector& Vector::operator=(const Vector &otherVector)
{
    assert(mSize == otherVector.mSize);
    for(int i=0; i<mSize; i++)
    {
        mData[i] = otherVector.mData[i];
    }
    return *this;
}

std::ostream& operator<<(std::ostream& output, Vector& v)
{
    output << std::endl;
    for(int i=0; i<v.Size(); i++)
    {
        std::cout << v[i] << std::endl;
    }
    return output; // return std::ostream so that we can have chain call of <<
}
