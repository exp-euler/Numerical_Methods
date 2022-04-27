# Runge-Kutta C++ Implementation

Runge-Kutta methods are numerical methods for solving initial value problems.
This repository provides a general implementation of such methods without
focusing in any one in particular. It includes a few methods already implemented
and one can very easily add more methods by only providing their Butcher
Tableaus. [This Wiki](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
article nicely presents the idea behind these methods,
how they are written in terms of a Butcher Tableau as well as provides the
tableaus of some famous RK methods.

# Usage

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
You should have successfully created the **`main`** file, which you can run in the
terminal to see the results:
```sh
user@user:~/Numerical_Methods/build$ ./main
```
