#!/bin/bash

#SBATCH --job-name=GL
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --time=60
#SBATCH --output=iks-%j.out

mkdir results_gl  ### All output goes here ###
mkdir results_gl2  ### All output goes here ###

pop=$1 ### population name ###
ref=$2 ### reference fasta file used in alignment ###


#	angsd -b ${pop}.bamlist -ref ${ref} -out results_gl/${pop} -uniqueOnly 2 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 -GL 2 -doGlf 4 -r 2

	angsd -b ${pop}.bamlist -ref ${ref} -out results_gl2/${pop} -uniqueOnly 2 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 -GL 2 -doGlf 1 -r 2

	angsd -glf results_gl2/${pop}.glf.gz -fai ${ref}.fai -nInd 10 -out results_gl2/${pop}.1 -doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1 -r 2
	angsd -glf results_gl2/${pop}.glf.gz -fai ${ref}.fai -nInd 10 -out results_gl2/${pop}.2 -doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1 -postCutoff 0.80 -r 2
	angsd -glf results_gl2/${pop}.glf.gz -fai ${ref}.fai -nInd 10 -out results_gl2/${pop}.3 -doMajorMinor 1 -doGeno 3 -doPost 1 -doMaf 1 -postCutoff 0.80 -r 2

