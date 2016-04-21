#!/bin/bash
## run with `sbatch -N[#machines] code/run.sbatch'
#SBATCH --time=8:0:0
#SBATCH --mem-per-cpu=2000
#SBATCH --output=results/%j.err
#SBATCH --constraint=ib
#SBATCH --exclusive

NOBS=$1
DESIGN=$2
SUPPORT=$3

RHO=$4
S2N=$5
DECAY=$6

export OUT=n$NOBS-$DESIGN-$SUPPORT/sim-rho$RHO-s2n$S2N-decay$DECAY   

LOG=results/$OUT-sbatch.log
echo `date` > $LOG
echo "job id $SLURM_JOBID" >> $LOG
echo "$SLURM_NNODES nodes" >> $LOG

module load parallel
npartask=$(($SLURM_NNODES * 16)) 
echo "$npartask tasks" >> $LOG

export para="parallel --delay .5 -j $npartask --joblog results/$OUT-para.log"
echo $para >> $LOG
export srun="srun --exclusive -N1 -c1 --ntasks=1"  
echo $srun >> $LOG
$para $srun "Rscript code/simrun.R {1} $NOBS $DESIGN $SUPPORT $RHO $S2N $DECAY >> results/$OUT-r.log" ::: {1..1000}
        
echo "done @" `date` >> $LOG



