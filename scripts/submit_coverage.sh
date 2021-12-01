#!/bin/sh

for S in 1 2
do
    qsub -v "SIMNUM=$S" coverage.pbs
done
