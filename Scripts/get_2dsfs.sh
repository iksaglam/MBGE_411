#!/bin/bash

#SBATCH --job-name=2dsfs
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --time=60
#SBATCH --output=iks-%j.out

pop1=$1
pop2=$2

                realSFS results_sfs/${pop1}.saf.idx results_sfs/${pop2}.saf.idx > results_sfs/${pop1}.${pop2}.2dsfs
                Rscript scripts/plot2DSFS.R results_sfs/${pop1}.${pop2}.2dsfs 10 10 $pop1 $pop2



