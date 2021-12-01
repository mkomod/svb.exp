#!/bin/sh

for S in 1 2 3 4 5 6
do
    qsub -v "SIMNUM=$S" simulations.pbs
done
