#!/bin/bash

#SBATCH --job-name=sfs
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --time=60
#SBATCH --output=iks-%j.out


mkdir results_sfs  ### All output goes here ###

pop=$1 ### population name ###
ref=$2 ### reference fasta file used in alignment ###
anc=$3 ### ancestral sequence file ### 


                angsd -b ${pop}.bamlist -ref $ref -anc $anc -out results_sfs/${pop} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -doCounts 1 -GL 2 -doSaf 1 -r 2
                realSFS results_sfs/${pop}.saf.idx > results_sfs/${pop}.sfs
                Rscript scripts/plotSFS.R results_sfs/${pop}.sfs $pop 0 results_sfs/${pop}.sfs.pdf



