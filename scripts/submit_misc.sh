#!/bin/sh

for S in 1 2 3 4 5
do
    qsub -v "DGP=$S" misc.pbs
done
