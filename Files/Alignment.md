
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

This will result in several indexes including the BWT index `~/course_content/week03_tutorial/references/isophya_contigs_CAYMY.fasta.bwt` which the alignments will be based on. The creation of the index only needs to be performed once and does not have to be recreated for every alignment job. 

### Aligning reads with BWA-MEM
Now that we have our indexes created, we can start aligning our reads (i.e. individual fastq files) to the reference genome. Let us first create a new directory where we will store our alignment files and cd (change directory) into it.

```Bash
mkdir alignments
cd alignments/
```

Next we will use the BWA-MEM algorithm to align one of our paired-end reads to the reference genome and output results as a SAM (Sequence Alignment/Map Format) file.

bwa mem ../references/isophya_contigs_CAYMY.fasta ../reads/CAYMY_002_R1.fastq.gz ../reads/CAYMY_002_R2.fastq.gz > CAYMY_002.sam

As before we can take a look at our SAM files using less.

```Bash
less -S CAYMY_002.sam
```

### SAM/BAM conversion, sorting and PCR clone (duplicate) markup/removal:

To save space we ideally want to transform our SAM file into a binary BAM format and sort by coordinates. Next we usually would like to work with only properly paired reads (i.e. where the forward and reversed reads mapped to the correct positions). Finally we want to remove/mark PCR duplicates using [picard-tools](https://broadinstitute.github.io/picard/), so they do not bias variant calling and genotyping and index our final BAM file.

```Bash
samtools view -bS CAYMY_002.sam > CAYMY_002.bam
samtools sort CAYMY_002.bam CAYMY_002_sorted.bam
```

To view our resulting bam files we can use the following command

```Bash
samtools view CAYMY_002_sorted.bam | less -S
```

Next we usually would like to work with only properly paired reads (i.e. where the forward and reversed reads mapped to the correct positions). Finally we want to remove/mark PCR duplicates using `picard-tools`, so they do not bias variant calling and genotyping and index our final BAM file.

```Bash
samtools view -b -f 0x2 CAYMY_002_sorted.bam > CAYMY_002_sorted_proper.bam
java -jar ~/bin/picard.jar MarkDuplicates INPUT=CAYMY_002_sorted_proper.bam OUTPUT=CAYMY_002_sorted_proper_rmdup.bam METRICS_FILE=CAYMY_002_metrics.txt VALIDATION_STRINGENCY=LENIENT  REMOVE_DUPLICATES=True
samtools index CAYMY_002_sorted_proper_rmdup.bam
```

For the purposes of this tutorial we are directly removing duplicates `REMOVE_DUPLICATES=True` from our resulting bam files (to primarily save space) but usually we only want to mark our duplicates and not remove them `REMOVE_DUPLICATES=False`.

### Alignment statistics

After aligning our reads to a reference genome a good first step would be to check the succes rate of our alignments. We can do this using a simple command in `samtools` called `flagstats`

```Bash
samtools flagstat A_CAMD04_sorted.bam



```



