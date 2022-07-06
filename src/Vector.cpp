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
