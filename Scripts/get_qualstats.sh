#!/bin/bash

#SBATCH --job-name=qs
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --time=60
#SBATCH --output=iks-%j.out


mkdir results_qs  ### All output goes here ###

infile=$1 ### list containing population names ###
ref=$2 ### reference fasta file used in alignment ###

for pop in `cat $infile`
do



        angsd -b ${pop}.bamlist -ref $ref -out results_qs/$pop -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 200 -minQ 0 -r 2
        Rscript scripts/plotQC.R results_qs/$pop

done
