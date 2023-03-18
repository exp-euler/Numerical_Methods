#include "TimeIntegration.hpp"
#include "Tableaus.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include <cstdlib>
#include <vector>
#include <iostream>
#include <mpi.h>
#include <time.h>

#include <omp.h>
#include <unistd.h>

int main(int argc, char*argv[])
{
    // Testing of the operations implemented in parallel with MPI

    int ProcNum;        // Number of available processes
    int ProcRank;       // Rank of current process

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    int m = 7000;
    int n = 3000;
    int k = 400;
    Matrix M(m,n);
    Matrix N(n,k);
    Matrix R(m,k);
    Matrix Rpar(m,k);

    double dtime;


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
                //M(i,j)=i+j;
                M(i,j)=1;
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
        dtime = omp_get_wtime();
        R = M.SerialMM(N);
        dtime = omp_get_wtime() - dtime;
        
        std::cout << "Result from serial multiplication: " << std::endl;
        std::cout << dtime << std::endl;
        //std::cout << R << std::endl;
    }

    // Wait for the serial multiplication to finish before continuing
    MPI_Barrier(MPI_COMM_WORLD);

    //Wpar = M*V;
    dtime = omp_get_wtime();
    Rpar = M*N;
    dtime = omp_get_wtime() - dtime;

    if(ProcRank == 0)
    {
        std::cout << "Result from parallel multiplication: " << std::endl;
        std::cout << dtime << std::endl;
        //std::cout << Rpar << std::endl;
    }


    MPI_Finalize();

    return 0;
}
