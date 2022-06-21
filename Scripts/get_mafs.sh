#!/bin/bash

#SBATCH --job-name=mafs
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --time=60
#SBATCH --output=iks-%j.out

mkdir results_pops_mafs  ### All output goes here ###

pop=$1 ### list containing population names ###
ref=$2 ### reference fasta file used in alignment ###

	angsd -b ${pop}.bamlist -ref ${ref} -out results_pops_mafs/${pop}.1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 -GL 1 -doGlf 1 -doMajorMinor 1 -doMaf 1 -r 2

	angsd -b ${pop}.bamlist -ref ${ref} -out results_pops_mafs/${pop}.2 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -doCounts 1 -GL 2 -doGlf 1 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 -r 2


