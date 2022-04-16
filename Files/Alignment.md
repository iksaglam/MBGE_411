
# Introduction to indexing and short read alignment


As a first step copy all content you will need in this exercise into your own directory and check that everything has been copied correcly.

```Bash
cd ~/my_directory/
cp -r ~/course_content/week03_tutorial/* ./
ls *
```

## Aligning reads to a reference genome

For aligning reads to a reference genoome we will use [BWA](http://bio-bwa.sourceforge.net/) (Burrows-Wheeler Aligner) which is based on the [Burrows-Wheeler Transformed](https://academic.oup.com/bioinformatics/article/25/14/1754/225615?login=true) index of the reference genome. It performs gapped alignment for single-end reads, supports paired-end mapping, generates mapping quality and can give multiple hits if called for.

To start let us create an empty shell script where we will store all our commands and later use to execute our pipeline.

```Bash
vim align_pe_reads.sh

#!/bin/bash

#SBATCH --job-name=align
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=mid
#SBATCH --mem=120G
#SBATCH --time=720
```

### Creating BWA-MEM index
Like all other alignment tools for genome wide short reads the first step is to index the reference genome. BWA indexes the genome with an FM Index based on the Burrows-Wheeler Transform enabling memory efficient mapping.


The basic command for indexing the genome using BWA which we will add to our shell script is:
```Bash
#!/bin/bash

#SBATCH --job-name=align
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=mid
#SBATCH --mem=120G
#SBATCH --time=720

ref=$1

bwa index -a bwtsw $ref
```
Now let us execute the script for one of the genomes using the following command:

```Bash
sbatch align_pe_reads.sh references/isophya_contigs_CAYMY.fasta
```

This will result in several indexes including the BWT index ~/course_content/week03_tutorial/references/isophya_contigs_CAYMY.fasta.bwt which the alignments will be based on. The creation of the index only needs to be performed once and does not have to be recreated for every alignment job. 



