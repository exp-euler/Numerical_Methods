#include "Matrix.hpp"
#include "Vector.hpp"
#include "MatrixFunctions.hpp"
#include <ctime>

#ifdef GPU
#include <driver_types.h>
__global__ void MM1 (float *a, float *b, float *c, int n)
{
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    float temp = 0;
    for(int k=0; k<n;k++)
        temp += a[i*n+k] * b[k*n+j];
    c[i*n+j] = temp;
}
#endif

int main(){
    srand(static_cast<unsigned>(time(0)));
    int size = 16*16;
    int rand_range = 2;

    // Objects for the scalar CPU multiplication
    Matrix<float> a(size,size);
    Matrix<float> b(size,size);
    Matrix<float> c(size,size);

    // Objects for the GPU multiplication
    float *A, *B, *C;           // host copies
    float *d_A, *d_B, *d_C;     // device copies
    size_t mem_size = size * size * sizeof(float);

    cudaMalloc((void **)&d_A, mem_size);
    cudaMalloc((void **)&d_B, mem_size);
    cudaMalloc((void **)&d_C, mem_size);

    A = (float *)malloc(mem_size);
    B = (float *)malloc(mem_size);
    C = (float *)malloc(mem_size);

    float numA;
    float numB;
    for(int i=0; i<size; i++) {
        for(int j=0; j<size; j++) {
            numA = 1 + (rand() % rand_range);
            numB = 2 + (rand() % rand_range);
            a(i,j) = numA;
            b(i,j) = numB;
            A[i*size+j] = numA;
            B[i*size+j] = numB;
        }
    }

    // Copy inputs to device
    cudaMemcpy(d_A, A, mem_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, mem_size, cudaMemcpyHostToDevice);

    // GPU multiplication
    int p=16; int q=16;
    dim3 grid (size/q, size/p);
    dim3 block (q,p);
    MM1 <<<grid, block>>> (d_A, d_B, d_C, size);

    cudaMemcpy(C, d_C, mem_size, cudaMemcpyDeviceToHost);

    // Serial CPU multiplication
    c = a*b;

    // Compare results
    for(int i=0; i<size; i++) {
        for(int j=0; j<size; j++) {
            if(C[i*size+j] - c(i,j) != 0) {
                std::cout << "Results differ" << std::endl;
                // Cleanup
                free(A); free(B); free(C);
                cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
            }
            return 0;
        }
    }


    // Cleanup
    free(A); free(B); free(C);
    cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);

    return 0;
}
