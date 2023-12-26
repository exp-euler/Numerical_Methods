#!/bin/bash

read -p "Initial number of steps: " init_steps
steps=$init_steps
read -p "Total number of runs: " num_runs
read -p "Order of the method used: " method_order
read -p "Executable to run: " executable
read -p "Problem solved: " problem

mkdir TEST_RUN && cd ./TEST_RUN
../../build/bin/$executable -N $steps
echo "Run with $steps steps."

for (( i=1; i<=$num_runs; i++ ))
do
    steps=$(( $steps * 2 ))
    ../../build/bin/$executable -N $steps
    echo "Run with $steps steps."
done

cd ..
python3 test_RK_order.py TEST_RUN $problem $method_order $init_steps $num_runs 1
