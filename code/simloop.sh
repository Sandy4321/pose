#!/bin/bash

for NOBS in 100 1000; do
	for DESIGN in "binary" "continuous"; do
		for SUPPORT in "sparse" "dense"; do
			mkdir results/n$NOBS-$DESIGN-$SUPPORT
			for RHO in 0 0.5 0.9; do
    			for S2N in 0.5 1 2; do
        			for DECAY in 10 50 100 200; do
        				echo $NOBS $DESIGN $SUPPORT $RHO $S2N $DECAY 
        				#bash code/run.sbatch $NOBS $DESIGN $SUPPORT $RHO $S2N $DECAY
            			sbatch -Jn$NOBS-$DESIGN-$SUPPORT-r$RHO-s$S2N-d$DECAY -N2 code/run.sbatch $NOBS $DESIGN $SUPPORT $RHO $S2N $DECAY
            		done
            	done
            done
        done
    done
done


