#!/bin/bash
## run with `sbatch -N[#machines] code/run.sbatch'
#SBATCH --time=12:0:0
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=2000
#SBATCH --output=%j.err
#SBATCH --constraint=ib
#SBATCH --exclusive

module load parallel

RHO=0.9
S2N=0.5
OUT=sim-rho$RHO-s2n$S2N
mkdir -p results/$OUT

LOG=results/$OUT/an_sbatch_$SLURM_JOBID.log
echo `date` > $LOG
echo "$SLURM_NNODES nodes" >> $LOG

export para="parallel --delay .2 -j $SLURM_NTASKS --joblog results/$OUT/para.log"
export srun="srun --exclusive -N1 -c1 --ntasks=1"  

$para $srun "Rscript code/simrun.R {1} $RHO $S2N $OUT > results/$OUT/r.log" ::: {1..1000}

echo "done @" `date` >> $LOG


