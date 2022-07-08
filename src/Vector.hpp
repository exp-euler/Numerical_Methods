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
        double& front();
        double& operator[](int i);
        Vector& operator=(const Vector& otherVector);
        friend std::ostream& operator<<(std::ostream& output, Vector& v);
};

#endif
