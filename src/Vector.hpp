#ifndef VECTORHEADERDEF
#define VECTORHEADERDEF

#include<vector>

class Vector
{
    private:
        std::vector<double> mData;
        int mSize;
    public:
        Vector(const Vector& otherVector);
        Vector(int size);
        //~Vector();
        int Size() const;
};

#endif
