# Numerical Solutions of Initial Value Problems

This repository aims to implement numerical time integrators of
Runge-Kutta type in a unified way. The main goal is to be able to perform
time integration of an ODE system in an efficient (at some point
parallelized on a GPU) manner that allows for automatic time-step
adjustment.

## Runge-Kutta C++ Implementation

Runge-Kutta (RK) methods are numerical methods for solving initial value problems.
A general implementation of such methods is provided, without
focusing in any one in particular. It includes a few methods already
implemented, but one can very easily add more methods by only providing
their Butcher
Tableaus. [This Wiki](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
article nicely presents the idea behind these methods,
how they are written in terms of a Butcher Tableau as well as provides the
tableaus of some famous RK methods.

## Exponential Runge-Kutta C++ Implementation

Exponential Runge-Kutta (ERK) methods are similar to the classical RK
methods, with the exception that they treat the linear part of the
ODE exactly. This is done in order to be able to filter out the
numerical stiffness of the ODE and allow for the use of a much larger
time-step in comparison to classical explicit RK methods.\

In my master thesis, I worked on exploring a better understanding of
numerical
stiffness (which can be found [here](https://arxiv.org/pdf/2305.12488.pdf)) as well as designed a robust
embedded ERK pair that allows for adjustable time-step (found [here](https://arxiv.org/pdf/2303.12139.pdf)).
The second aim of this repository is to have a C++ implementation of
the works mentioned above, which would enable the
testing of new ideas on the subject of numerical stiffness.\

One of the difficulties of implementing ERK methods is the accurate and
fast calculation of the so called \phi functions. The norm in the
literature is to use a form of the scaling and squaring algorithm coupled
with a Pade approximation method, but I am experimenting with
Taylor approximations instead.


## Usage

To run the project, first create a **`build/`** directory by running:
```sh
user@user:~/Numerical_Methods$ mkdir build
user@user:~/Numerical_Methods$ cd build
```
and then run cmake like below:
```sh
user@user:~/Numerical_Methods/build$ cmake ../
user@user:~/Numerical_Methods/build$ make -j
```
You should have successfully created **`main`** files, which you can run
in the terminal to see the results:
```sh
user@user:~/Numerical_Methods/build$ ./main
```

As described above, no third party library needs to be installed in
order for the compilation and run to succeed. In addition, the code can
be compiled with OpenMP , MPI or CUDA support
(by altering **`CMAKE_BUILD_TYPE`**). As this is a numerical
algorithm, I have also implemented a Matrix class and a Vector class.
Nevertheless, it is very straightforward and up to the user to
make a few simple changes in a few of the files in order to use
the Eigen package instead.
