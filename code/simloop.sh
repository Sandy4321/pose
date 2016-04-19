#!/bin/bash

for RHO in 0 0.5 0.9; do
    for S2N in 0.5 1 2; do
        for DECAY in 10 50 100 200; do
            sbatch -Jr$RHO-s$S2N-d$DECAY -N2 code/run.sbatch $RHO $S2N $DECAY
        done
    done
done


