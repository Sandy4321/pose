#!/bin/bash
## run with `sbatch -N[#machines] code/run.sbatch'
#SBATCH --time=12:0:0
#SBATCH --mem-per-cpu=2000
#SBATCH --output=%j.err
#SBATCH --constraint=ib
#SBATCH --exclusive

module load parallel

npartask=$(($nn * 16)) 

for RHO in 0 0.5 0.9; do
    for S2N in 0.5 1 2; do
        for DECAY in 10 50 100 200; do

            OUT=sim-rho$RHO-s2n$S2N-decay$DECAY        
            LOG=results/$OUT-sbatch.log
            echo `date` > $LOG
            echo "$SLURM_NNODES nodes" >> $LOG
            echo "$npartask tasks" >> $LOG

            export para="parallel --delay .2 -j $npartask --joblog results/$OUT-para.log"
            echo $para >> $LOG
            export srun="srun --exclusive -N1 -c1 --ntasks=1"  
            echo $srun >> $LOG
        
            $para $srun "Rscript code/simrun.R {1} $RHO $S2N $DECAY $OUT > results/$OUT-r.log" ::: {1..1000}
        
            echo "done @" `date` >> $LOG
        done
    done
done


