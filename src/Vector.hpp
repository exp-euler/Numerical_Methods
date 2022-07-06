#ifndef VECTORHEADERDEF
#define VECTORHEADERDEF

#include<vector>
#include<iostream>

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
        double& operator[](int i);
        friend std::ostream& operator<<(std::ostream& output, Vector& v);
};

#endif
