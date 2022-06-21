#!/bin/bash

#SBATCH --job-name=geno
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --time=60
#SBATCH --output=iks-%j.out

mkdir results_pops_geno  ### All output goes here ###

infile=$1 ### list containing population names ###
ref=$2 ### reference fasta file used in alignment ###

for pop in `cat $infile`
do

echo "#!/bin/bash 

#SBATCH --job-name=geno
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=short
#SBATCH --time=60
#SBATCH --output=iks-%j.out


	angsd -b ${pop}.bamlist -ref ${ref} -out results_pops_geno/${pop} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 2 -doCounts 1 -GL 2 -doGlf 1 -doMajorMinor 1 -doGeno 5 -doPost 1 -doMaf 1 -postCutoff 0.80 -r 2 -sites snp.txt" > ${pop}.sh

	sbatch ${pop}.sh

done

