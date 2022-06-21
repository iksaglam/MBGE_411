#!/bin/bash

#SBATCH --job-name=PBS
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --time=60
#SBATCH --output=iks-%j.out

mkdir results_fst

pop1=$1
pop2=$2
pop3=$3

### Calculate index for fst and pbs analysis###

	realSFS fst index results_sfs/${pop1}.saf.idx results_sfs/${pop2}.saf.idx results_sfs/${pop3}.saf.idx -sfs results_sfs/${pop1}.${pop2}.2dsfs -sfs results_sfs/${pop1}.${pop3}.2dsfs -sfs results_sfs/${pop2}.${pop3}.2dsfs -whichFST 1 -fstout results_fst/CEU.pbs

### Calculate global fst and pbs statistics ###

	realSFS fst stats2 results_fst/CEU.pbs.fst.idx -win 50000 -step 10000 -whichFST 1 > results_fst/CEU.pbs.fst.txt
	Rscript scripts/plotPBS.R results_fst/CEU.pbs.fst.txt results_fst/CEU.pbs.fst.txt.pdf 135 139

