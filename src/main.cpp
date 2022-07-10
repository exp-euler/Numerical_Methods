#include "TimeIntegration.hpp"
#include "Tableaus.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include <cstdlib>
#include <vector>
#include <iostream>
#include <mpich/mpi.h>
#include <time.h>

#include <omp.h>
#include <unistd.h>

// Demonstration of specifying a Runge-Kutta method to solve an equation like:
//
//                     dy/dt = F(t,y),    y(0) = y0
//
// where the function F(t,y) can be arbitrary and is referred as Right Hand Side.
double RHS(double t, double y)
{
    return 1+t;
}

int main(int argc, char*argv[])
{
    /*
    // Testing of the time integration part of the code (serial for now)
    TimeIntegration equation;
    equation.Solve(ClassicalRK::Euler(), &RHS, (0.0+1.0)/100, 2.0, 0, 1, 100);
    const std::vector<double> &sol = equation.getY();
    */

    // Testing of the operations implemented in parallel with MPI

    int ProcNum;        // Number of available processes
    int ProcRank;       // Rank of current process

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    int m = 7;
    int n = 3;
    int k = 4;
    Matrix M(m,n);
    Matrix N(n,k);
    Matrix R(m,k);
    Matrix Rpar(m,k);


    Vector V(n);
    Vector W(m);
    Vector Wpar(m);

    // Use time as seed for random number generator
    srand(time(0));

    if(ProcRank == 0)
    {
        for(int i=0; i<m; i++)
        {
            for(int j=0; j<n; j++)
            {
                M(i,j)=i*n+j;
                //M(i,j)=rand() % 100;
            }
        }

        for(int i=0; i<n; i++)
        {
            for(int j=0; j<k; j++)
            {
                N(i,j)=i*n+j;
                //M(i,j)=rand() % 100;
            }
        }

        for(int i=0; i<n; i++)
        {
            //V[i] = i;
            V[i] = rand() % 100;
        }
    }

    // Wait for the matrix and vector to be initialized before continuing
    MPI_Barrier(MPI_COMM_WORLD);

    if(ProcRank == 0)
    {
        //W = M.SerialMV(V);
        R = M.SerialMM(N);
        
        std::cout << "Result from serial multiplication: " << std::endl;
        std::cout << R << std::endl;
    }

    // Wait for the serial multiplication to finish before continuing
    MPI_Barrier(MPI_COMM_WORLD);

    //Wpar = M*V;
    Rpar = M*N;

    if(ProcRank == 1)
    {
        std::cout << "Result from parallel multiplication: " << std::endl;
        std::cout << Rpar << std::endl;
    }


    MPI_Finalize();

    /*
    #pragma omp parallel
    {
        usleep(5000 * omp_get_thread_num());
        std::cout << "Hellow World... from thread = " << omp_get_thread_num();
        std::cout << std::endl;
    }
    */

    return 0;
}
