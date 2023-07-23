#!/bin/bash

read -p "Initial number of steps: " steps
read -p "Total number of runs: " num_runs
../../../build/bin/vectorIVP -N $steps
echo "Run with $steps steps."

for (( i=1; i<=$num_runs; i++ ))
do
    steps=$(( $steps * 2 ))
    ../../../build/bin/vectorIVP -N $steps
    echo "Run with $steps steps."
done
