# An Introduction to Comparing Genomes: Selection for lactose tolerance in Europeans 


### Objectives:
- Are functional alleles causing lactose tolerance found in higher frequencies in European populations?
- Can we identify signatures of natural selection around the LCT (Lactase-phlorizin hydrolase) gene in European populations?

### Action Plan:
- Estimate allele frequencies for known variants in African (LWK, YRI) and European (CEU) populations
- Perform a sliding windows scan based on allele frequency differentiation


All data, scripts and results in this tutorial can be found in the following directory: `/kuacc/users/mbge411/hpc_run/2024SpringMBGE411/week04_tutorial`


## Practical 1: Estimating allele frequencies and variants in African (LWK, YRI) and European (CEU) populations

### In this practical we will learn how to do:
- Call genotypes
- Estimate allele frequencies
- Determine variant/polymorphic sites

To accomplish these goals, we will use the program [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) (Analysis of Next Generation Sequencing Data) developed by [Thorfinn Korneliussen](https://scholar.google.com/citations?user=-YNWF4AAAAAJ&hl=en) and [Anders Albrechtsen](https://scholar.google.dk/citations?user=20oVxFsAAAAJ&hl=en) (and many other contributors) at the University of Copenhagen and University of California Berkeley. More information about the program and implemented methods can be found [here](https://www.ncbi.nlm.nih.gov/pubmed/25420514).


### Basic Filtering with ANGSD

First let’s take a look at our options in ANGSD:

```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd
```

We should see something like this:

```Bash
...
Overview of methods:
        -GL             Estimate genotype likelihoods
        -doCounts       Calculate various counts statistics
        -doAsso         Perform association study
        -doMaf          Estimate allele frequencies
        -doError        Estimate the type specific error rates
        -doAncError     Estimate the errorrate based on perfect fastas
        -HWE_pval               Est inbreedning per site or use as filter
        -doGeno         Call genotypes
        -doFasta        Generate a fasta for a BAM file
        -doAbbababa     Perform an ABBA-BABA test
        -sites          Analyse specific sites (can force major/minor)
        -doSaf          Estimate the SFS and/or neutrality tests genotype calling
        -doHetPlas      Estimate hetplasmy by calculating a pooled haploid frequency
        Below are options that can be usefull
        -bam            Options relating to bam reading
        -doMajorMinor   Infer the major/minor using different approaches
        -ref/-anc       Read reference or ancestral genome
        -doSNPstat      Calculate various SNPstat
        -cigstat        Printout CIGAR stat across readlength
        many others
For information of specific options type:
        ./angsd METHODNAME eg
                ./angsd -GL
                ./angsd -doMaf
                ./angsd -doAsso etc
                ./angsd sites for information about indexing -sites files
Examples:
        Estimate MAF for bam files in 'list'
                './angsd -bam list -GL 2 -doMaf 2 -out RES -doMajorMinor 1'
```
If the input file is in BAM format, the possible options are:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd -bam
...
parseArgs_bambi.cpp: bam reader:
        -bam/-b         (null)  (list of BAM/CRAM files)
        -i              (null)  (Single BAM/CRAM file)
        -r              (null)  Supply a single region in commandline (see examples below)
        -rf             (null)  Supply multiple regions in a file (see examples below)
        -remove_bads    1       Discard 'bad' reads, (flag >=256)
        -uniqueOnly     0       Discards reads that doesnt map uniquely
        -show           0       Mimic 'samtools mpileup' also supply -ref fasta for printing reference column
        -minMapQ        0       Discard reads with mapping quality below
        -minQ           13      Discard bases with base quality below
        -trim           0       Number of based to discard at both ends of the reads
        -trim           0       Number of based to discard at 5 ends of the reads
        -trim           0       Number of based to discard at 3 ends of the reads
        -only_proper_pairs 1    Only use reads where the mate could be mapped
        -C              0       adjust mapQ for excessive mismatches (as SAMtools), supply -ref
        -baq            0       adjust qscores around indels (as SAMtools), supply -ref
        -checkBamHeaders 1      Exit if difference in BAM headers
        -doCheck        1       Keep going even if datafile is not suffixed with .bam/.cram
        -downSample     0.000000        Downsample to the fraction of original data
        -nReads         50      Number of reads to pop from each BAM/CRAMs
        -minChunkSize   250     Minimum size of chunk sent to analyses
Examples for region specification:
                chr:            Use entire chromosome: chr
                chr:start-      Use region from start to end of chr
                chr:-stop       Use region from beginning of chromosome: chr to stop
                chr:start-stop  Use region from start to stop from chromosome: chr
                chr:site        Use single site on chromosome: chr
```
Now let us define filtering options
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd -b human_pop.bamlist -ref references/hum_ch2_ref.fa -out results_qs/human_pop -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 200 -minQ 0 -r 2
```
These filters will retain only uniquely mapping reads, not tagged as bad, considering only proper pairs, without trimming, and adjusting for indel/mapping (as in samtools). 
`-C 50` reduces the effect of reads with excessive mismatches, while `-baq 1` computes base alignment quality (BAQ) as explained [here](http://samtools.sourceforge.net/mpileup.shtml) and is used to rule out false SNPs close to INDELS.
We can take a look out our results:
```Bash
ls results/*
...
results/human_pop.arg  Results/human_pop.depthGlobal results/human_pop.depthSample results/human_pop.qs
...
```
Let us compute the percentiles of these distributions (and visualize them) in order to make informative decisions on threshold values for filtering. We can use an Rscript to do this
```Bash
Rscript scripts/plotQC.R results_qs/human_pop.qs
```
**QUESTION:** Which values would you choose as sensible thresholds on quality score and depth (minimum and maximum)?
It is also always good practice to remove sites where a fraction (usually half) of the individuals have no data. This is achieved by the `-minInd` option
So a good command line for filtering data could look like this:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd -b human_pop.bamlist -ref references/hum_ch2_ref.fa -out Results/human_pop -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 6 -setMaxDepth 30 -doCounts 1 -r 2
```
Now we are ready to calculate the genotype likelihoods at each site for each individual.
To do this we need to specify which genotype likelihood model to use:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd -GL
...
-GL=0:
        1: SAMtools
        2: GATK
        3: SOAPsnp
        4: SYK
        5: phys
        6: Super simple sample an allele type GL. (1.0,0.5,0.0)
        -trim           0               (zero means no trimming)
        -tmpdir         angsd_tmpdir/   (used by SOAPsnp)
        -errors         (null)          (used by SYK)
        -minInd         0               (0 indicates no filtering)
Filedumping:
        -doGlf  0
        1: binary glf (10 log likes)    .glf.gz
        2: beagle likelihood file       .beagle.gz
        3: binary 3 times likelihood    .glf.gz
        4: text version (10 log likes)  .glf.gz
```
As we can see there are different models we can use. For now we will choose the first one: SAMtools
Let us estimate genotype likelihoods (GL) for the European populations. A possible command line to estimate GLs might be:
```Bash
/userfiles/iksaglam/bin/angsd/angsd -b CEU.bamlist -ref references/hum_ch2_ref.fa -out results_GL/CEU -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 -GL 1 -doGlf 4 -r 2
```
where we specify:
- GL 1: genotype likelihood model as in SAMtools
- doGlf 4: output in text format
- doCounts 1: Count the number A,C,G,T. All sites, All samples
**QUESTION:** What are the output files? What's the information inside them?
```Bash
ls results_gl/CEU*
less -S results_gl/CEU.arg
zless -S results_gl/CEU.glf.gz
```
First 2 columns of the output are the genomic positions and the remaining lines are the 10 genotype likelihoods for each individual in the usual ordering: `AA AC AG AT CC CG CT GG GT TT`.
We can quickly count the number of columns in our outfile to confirm this:
```Bash
zless results_qs/CEU.glf.gz | awk '{print NF}' | sort -nu | tail -n 1
102
```
### Genotype calling:
We are now ready to call genotypes from sequencing data. We calculate genotype probabilities at each site for each individual using the genotype likelihoods calculated above.
In ANGSD, the option to call genotypes is `-doGeno`:
```Bash
/userfiles/iksaglam/bin/angsd/angsd -doGeno
...
-doGeno 0
        1: write major and minor
        2: write the called genotype encoded as -1,0,1,2, -1=not called
        4: write the called genotype directly: eg AA,AC etc
        8: write the posterior probability of all possible genotypes
        16: write the posterior probability of called genotype
        32: write the posterior probabilities of the 3 genotypes as binary
        -> A combination of the above can be chosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3
        -postCutoff=0.333333 (Only genotype to missing if below this threshold)
        -geno_minDepth=-1       (-1 indicates no cutoff)
        -geno_maxDepth=-1       (-1 indicates no cutoff)
        -geno_minMM=-1.000000   (minimum fraction af major-minor bases)
        -minInd=0  (only keep sites if you call genotypes from this number of individuals)
        NB When writing the posterior the -postCutoff is not used
        NB geno_minDepth requires -doCounts
        NB geno_maxDepth requires -doCounts
```
Let us set `-doGeno 4`, so our genotypes are coded as AA,AC etc. We also want to print the major and minor alleles so we set `-doGeno 5`.
To calculate the posterior probability of genotypes we need to define a model.
```Bash
/userfiles/iksaglam/bin/angsd/angsd -doPost
...
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
        3: Using SFS as prior (still in development)
        4: Using reference panel as prior (still in development), requires a site file with chr pos major minor af ac an
...
```
For now let us choose `-doPost 2`, meaning we have chosen a uniform prior.
Furthermore, to calculate genotype posterior probabilities we need to assign the major and minor alleles (if biallelic).
```Bash
/userfiles/iksaglam/bin/angsd/angsd -doMajorMinor
...
        -doMajorMinor   0
        1: Infer major and minor from GL
        2: Infer major and minor from allele counts
        3: use major and minor from a file (requires -sites file.txt)
        4: Use reference allele as major (requires -ref)
        5: Use ancestral allele as major (requires -anc)
        -rmTrans: remove transitions 0
        -skipTriallelic 0
```
A typical command for genotype calling is (assuming we are analyzing CEU samples):
```Bash
/userfiles/iksaglam/bin/angsd/angsd -b CEU.bamlist -ref references/hum_ch2_ref.fa -out results_gl/CEU -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 -GL 2 -doGlf 1
/userfiles/iksaglam/bin/angsd/angsd -glf results_gl/CEU.glf.gz -fai references/hum_ch2_ref.fa.fai -nInd 10 -out results_gl/CEU  -doMajorMinor 1 -doGeno 5 -doPost 2 -doMaf 1
```
Let us take a look at the output file:
```Bash
less -S results_gl/CEU.geno.gz
```
The columns are: chromosome, position, major allele, minor allele and genotypes (AA,AC etc.).
**QUESTION:** How many sites have at least one missing genotype?
```Bash
zless results_gl/CEU.geno.gz | grep NN - | wc -l
400571
```
Why is that?
You can control how to set a site as “missing genotype” in an individual when their confidence is low with `-postCutoff`. For instance, we can set a site as “missing genotype” when the highest genotype posterior probability at that site for that individual is below `0.95`:
```Bash
/userfiles/iksaglam/bin/angsd/angsd -glf results_gl/CEU.glf.gz -fai references/hum_ch2_ref.fa.fai -nInd 10 -out results_gl/CEU  -doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1 -postCutoff 0.95
```
How many sites do we have in total? How many sites have at least one missing genotype now?
```Bash
zless results_gl/CEU.geno.gz | wc -l
300586
zless results_gl/CEU.geno.gz | grep NN - | wc -l
300586
```
Let us now use an informative prior `-doPost 1` to  calculate genotype posterior probabilities by assuming HWE and using frequency of observed alleles. The command line for this would be:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd -glf results_gl/CEU.glf.gz -fai references/hum_ch2_ref.fa.fai -nInd 10 -out results_gl/CEU_2 -doMajorMinor 1 -doGeno 5 -doPost 1 -doMaf 1
```
When we look at the results we can see that all sites have been called and we no longer have any missing sites.
```Bash
zless results_gl/CEU_2.geno.gz | grep NN - | wc -l
0
```
### Objective 1: Are functional alleles causing lactose tolerance found in higher frequencies in European populations?
We now have enough knowledge to build an analytical pipeline for answering our first question: Are functional alleles for lactose tolerance found at higher frequencies in European populations. 
Our SNPs of interest are located in the MCM6 gene which lies upstream of the LCT gene and contains two regulatory regions for LCT within its introns. SNP designation and positions of these variants can be found in `functional_snps.txt`. 
Let us look at this file and breakdown the property of these SNPs:
```Bash
cat functional_snps.txt
SNP	Chr	Pos
rs4988235	2	136608646
rs182549	2	136616754
```
**rs4988235:** rs4988235 (G/A) is associated with hypolactasia, more commonly known as lactose intolerance in European Caucasian populations. In these populations, the rs4988235(A) allele is both the more common allele and the one associated with lactase persistence; individuals who are rs4988235(G;G) are likely to be lactose intolerant.
(G;G = likely to be lactose intolerant as an adult) (A;G = likely to be able to digest milk as an adult) (A;A = can digest milk).
**rs182549:** rs182549 (C/T) is associated with hypolactasia, more commonly known as lactose intolerance in European Caucasian populations. In these populations, the rs182549(C) allele is both the more common allele and the one associated with lactose intolerance. (C;C = possibly lactose intolerant) (C;T = Can digest milk) (T;T = Can digest milk).
In ANGSD we can restrict our analyses on a subset of positions of interest using the `-sites` option. The file with these positions need to be formatted as (chromosome positions).
```Bash
sed 1d  functional_snps.txt | awk '{print $2, $3}' > snp.txt
```
We need to index this file in order for ANGSD to process it.
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd sites index snp.txt
```
We will now determine genotypes in each population and then calculate allele frequencies based on called genotypes. We can do this with the following code:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd -b ${pop}.bamlist -ref ${ref} -out results_pops_geno/${pop} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 2 -doCounts 1 -GL 2 -doGlf 1 -doMajorMinor 1 -doGeno 5 -dovcf 1 -doPlink 2 -doPost 1 -doMaf 1 -postCutoff 0.80 -r 2 -sites snp.txt
```
In the above code we do the following:
- use the SAMtools genotype likelihood model `-GL 1`
- calculate genotype posterior probabilities using a HWE-based prior `-doPost 1`
- filter out bases with a quality score less than 20 `-minQ 20`
- filter our reads with a mapping quality score less than 20 `-minMapQ 20`
- use ony sites where you have at least two samples with data `-mindInd 2`
- do not set any filtering based on min and max depth
- use `-doMajorMinor 1` and `-doMaf 1` options
- set genotypes missing if genotype probability is less than 0.80 `-postCutoff 0.80`
- use option `-sites snp.txt` to restrict the analysis only on selected sites
When the analysis is complete, open the output files and calculate allele frequency by counting genotypes. Can you comment on these results? Do you see any allele frequency differentiation between European and African populations?
### Estimation of allele frequencies and SNP calling
We now want to estimate allele frequencies at each site. In other words, at each site we want to estimate (or count) how many copies of different alleles (two in case of biallelic variants) we observe in our sample (across all sequenced individuals).
ANGSD has an option to estimate allele frequencies taking into account data uncertainty from genotype likelihoods:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd -doMaf
...
-doMaf  0 (Calculate persite frequencies '.mafs.gz')
        1: Frequency (fixed major and minor)
        2: Frequency (fixed major unknown minor)
        4: Frequency from genotype probabilities
        8: AlleleCounts based method (known major minor)
        NB. Filedumping is supressed if value is negative
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
        3: Using SFS as prior (still in development)
        4: Using reference panel as prior (still in development)
Filters:
        -minMaf         -1.000000       (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
        -rmTriallelic   0.000000        (Remove sites with a pvalue lower)
Extras:
        -ref    (null)  (Filename for fasta reference)
        -anc    (null)  (Filename for fasta ancestral)
        -eps    0.001000 [Only used for -doMaf &8]
        -beagleProb     0 (Dump beagle style postprobs)
        -indFname       (null) (file containing individual inbreedcoeficients)
        -underFlowProtect       0 (file containing individual inbreedcoeficients)
NB These frequency estimators requires major/minor -doMajorMinor
```
Since we will be counting alleles the estimation of allele frequencies requires the specification of how to assign the major and minor alleles (if biallelic).
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd -doMajorMinor
...
        -doMajorMinor   0
        1: Infer major and minor from GL
        2: Infer major and minor from allele counts
        3: use major and minor from a file (requires -sites file.txt)
        4: Use reference allele as major (requires -ref)
        5: Use ancestral allele as major (requires -anc)
        -rmTrans: remove transitions 0
        -skipTriallelic 0
```
A command line to estimate allele frequencies for populations might look like this:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd -b ${pop}.bamlist -ref ${ref} -out results_pops_mafs/${pop} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1
```
This analysis will return an output file ending in `.mafs.gz`. Let us take a look:
```Bash
zless results_pops_mafs/CEU.mafs.gz | head
chromo	position	major	minor	ref	knownEM	nInd
2	119999922	C	A	C	0.000002	5
2	119999923	A	C	A	0.000002	5
2	119999924	G	A	G	0.000001	5
2	119999925	G	A	G	0.000002	6
2	119999926	C	A	C	0.000003	7
2	119999927	C	A	C	0.000003	7
2	119999928	T	A	T	0.000003	7
2	119999929	G	A	G	0.000004	8
2	119999930	G	A	G	0.000004	8
```
The columns are: `chromosome`, `position`, `major allele`, `minor allele`, `reference allele`, `minor allele frequency (knownEM)` and the `number of individuals` with data. Notice that many sites have a very low minor allele frequency (nearly zero). This indicates that the site is monomorphic. We may be interested in looking at allele frequencies only for sites that are actually variable in our sample. Therefore we want to perform a SNP calling.
There are two main ways to call SNPs using ANGSDs:
```Bash
-minMaf         0.000000        (Remove sites with MAF below)
-SNP_pval       1.000000        (Remove sites with a pvalue larger)
```
We can filter our allele frequency file so that only sites with a minor allele frequency over 0.05 are kept `-minMaf 0.05` or/and keep only those sites whose probability of being variable is over a specified p value `-SNP_pval 1e-12`.
Let us now do SNP calling using the following code for our populations:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd -b ${pop}.bamlist -ref ${ref} -out results_pops_mafs/${pop} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 -GL 2 -doGlf 1 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -SNP_pval 1e-12 -r 2
```
What are the minor allele frequencies for the functional SNPs `rs4988235 and rs182549` controlling lactose tolerance? Are they similar to the ones we calculated using called genotypes? Can we find these SNPs in the African populations?
## Identifying signatures of natural selection around the LCT (Lactase-phlorizin hydrolase) gene in European populations.
Here we will perform a scan for positive selection by calculating PBS (population branch statistic) in “genomic windows” from low-depth data. Specifically, we will estimate:
- site frequency spectrum of populations
- genetic differentiation between populations
### Allele frequency differentiation
The joint Site Frequency Spectrum (SFS) can be considered as a summary statistics for the level of genetic differentiation between populations. Here we will learn how to calculate the `SFS` of populations and `2D-SFS` between populations using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD). We will also see how we can use calculated SFSs as prior information to estimate PBS and other summary statistics from genome wide data. 
Our final goal is to detect signatures of selection in our data, around the LCT gene in European populations. To compute `FST/PBS` we first need to estimate the marginal and joint sites frequency spectra for our populations.
One of the most important aspects of data analysis for population genetics is the estimate of the `Site Frequency Spectrum (SFS)`. `SFS` records the proportion of sites with different number of counts for the derived allele. It can be folded or unfolded, and the latter case implies the use of an outgroup species to define the ancestral state. SFS is informative on the demography of the population or on selective events (when estimated at a local scale).
To calculate `SFS` of populations using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) we will first compute genotype likelihoods. Then from these quantities [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) computes posterior probabilities of `Sample Allele Frequency (SAF)`, for each site. Finally, an estimate of the `SFS` is computed.
These steps can be accomplished in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) using `-doSaf 1/2` options and the program `realSFS`.
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd -doSaf
...
-doSaf		0
	1: perform multisample GL estimation
	2: use an inbreeding version
	3: calculate genotype probabilities (use -doPost 3 instead)
	4: Assume genotype posteriors as input (still beta) 
	-doThetas		0 (calculate thetas)
	-underFlowProtect	0
	-fold			0 (deprecated)
	-anc			(null) (ancestral fasta)
	-noTrans		0 (remove transitions)
	-pest			(null) (prior SFS)
	-isHap			0 (is haploid beta!)
	-doPost			0 (doPost 3,used for accessing saf based variables)
NB:
	  If -pest is supplied in addition to -doSaf then the output will be posterior probability of the sample allele frequency for each site
```
We will compute the `SFS` for each population separately. We want to estimate the `unfolded SFS` so we will use our `reference genome` to determine `ancestral` and `derived` states. The following code can be used to calculate `SFS` of populations:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/angsd -b ${pop}.bamlist -ref $ref -anc $anc -out results_sfs/${pop} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -doCounts 1 -GL 2 -doSaf 1
```
Let us take a look at the output file.
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/misc/realSFS print results_sfs/CEU.saf.idx | less -S
```
These values represent the `sample allele frequency likelihoods` at each site. The first value (after the chromosome and position columns) is the likelihood of having 0 copies of the derived allele, the second indicates the probability of having 1 copy and so on. Note that these values are in log format and scaled so that the maximum is 0.
The next step would be to use these likelihoods and estimate the overall SFS. This is achieved by the program realSFS.
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/misc/realSFS
-> ---./realSFS------
	-> EXAMPLES FOR ESTIMATING THE (MULTI) SFS:
	-> Estimate the SFS for entire genome??
	-> ./realSFS afile.saf.idx 
	-> 1) Estimate the SFS for entire chromosome 22 ??
	-> ./realSFS afile.saf.idx -r chr22 
	-> 2) Estimate the 2d-SFS for entire chromosome 22 ??
	-> ./realSFS afile1.saf.idx  afile2.saf.idx -r chr22 
	-> 3) Estimate the SFS for the first 500megabases (this will span multiple chromosomes) ??
	-> ./realSFS afile.saf.idx -nSites 500000000 
	-> 4) Estimate the SFS around a gene ??
	-> ./realSFS afile.saf.idx -r chr2:135000000-140000000 
	-> Other options [-P nthreads -tole tolerence_for_breaking_EM -maxIter max_nr_iterations -bootstrap number_of_replications]
	-> See realSFS print for possible print options
	-> Use realSFS print_header for printing the header
	-> Use realSFS cat for concatenating saf files
	->------------------
	-> NB: Output is now counts of sites instead of log probs!!
	-> NB: You can print data with ./realSFS print afile.saf.idx !!
	-> NB: Higher order SFSs can be estimated by simply supplying multiple .saf.idx files!!
	-> NB: Program uses accelerated EM, to use standard EM supply -m 0
```
We can use the following command to estimate the SFS for each population:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/misc/realSFS results_sfs/${pop}.saf.idx > results_sfs/${pop}.sfs
```
Let us have a look at the output of one of the populations (CEU):
```Bash
cat results_sfs/CEU.sfs
```
The first value represents the expected number of sites with derived allele frequency equal to 0, the second column the expected number of sites with frequency equal to 1 and so on.
We can plot and visualize the SFS of our populations using a simple Rscript and the code below:
```Bash
Rscript scripts/plotSFS.R results_sfs/${pop}.sfs $pop 0 results_sfs/${pop}.sfs.pdf
```
What can you comment about the shape of the `SFS` for each population. How do they differ? How would you expect the shape of the SFS for European populations to look like? Do the results match your expectations?
Now, we need to estimate a `multi-dimensional SFS` between our populations. The `joint-SFS` between populations can be used for making inferences about the divergence of populations (time, migration rate and so on). However, here we are interested in estimating the `2D-SFS` as prior information for our `FST/PBS`.
An important issue when doing this is to be sure that we are comparing the exact same sites between populations. [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) does that automatically and considers only a set of overlapping sites.
We are performing `PBS` assuming `CEU (Europeans)` are the target population (i.e. target of selection), and `LWK` and `YRI` as reference populations. `2D-SFS` between populations can be computed using the following code:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/misc/realSFS results_sfs/LWK.saf.idx results_sfs/YRI.saf.idx > results_sfs/LWK.YRI.2dsfs
/kuacc/users/mbge411/hpc_run/bin/angsd/misc/realSFS results_sfs/LWK.saf.idx results_sfs/CEU.saf.idx > results_sfs/LWK.CEU.2dsfs
/kuacc/users/mbge411/hpc_run/bin/angsd/misc/realSFS results_sfs/YRI.saf.idx results_sfs/CEU.saf.idx > results_sfs/YRI.CEU.2dsfs
```
The output file is a flatten matrix, where each value is the count of sites with the corresponding joint frequency ordered as \[0,0\] \[0,1\] and so on.
```Bash
less -S results_sfs/LWK.CEU.sfs
```
Let us plot and look at our 2D-SFSs:
```Bash
Rscript scripts/plot2DSFS.R results_sfs/LWK.YRI.2dsfs LWK-YRI 10-10
Rscript scripts/plot2DSFS.R results_sfs/LWK.CEU.2dsfs LWK-CEU 10-10
Rscript scripts/plot2DSFS.R results_sfs/YRI.CEU.2dsfs YRI-CEU 10-10
```
Look at the `2D-SFS`, can we spot candidate regions with high differentiation, indicative of selection?
### Population Branch Statitics
Now that we have `2D-SFS` between all populations we can calculate a`llele frequency differentiation` using the `PBS (population branch statistic)` metric. Also, `CEU` is our target population, while `LWK` and `YRI` are reference populations. The `2D-SFS` of populations will be used as prior information for the `joint allele frequency probabilities` at each site. From these probabilities we will calculate the `population branch statistic (PBS)` using the `CEU` as `target` population and `LWK` and `YRI` as `reference` populations. Our goal is to `detect selection in CEU` in terms of allele frequency differentiation.
Specifically, we are computing a `sliding windows scan`, with `windows of 50kbp` and a `step of 10kbp`. This can be achieved using the following commands:
The first command will compute per-site FST indexes (please note the order of files):
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/misc/realSFS fst index results_sfs/LWK.saf.idx results_sfs/YRI.saf.idx results_sfs/CEU.saf.idx -sfs results_sfs/LWK.YRI.2dsfs -sfs results_sfs/LWK.CEU.2dsfs -sfs results_sfs/YRI.CEU.2dsfs -whichFST 1 -fstout results_fst/CEU.pbs
```
we can take a look at these values using the following command:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/misc/realSFS fst print results_fst/CEU.pbs.fst.idx | less -S
```
where columns are: `chromosome`, `position`, `(a), (a+b) values` for the three FST comparisons, where `FST is defined as a/(a+b)`. Note that FST on multiple SNPs is calculated as `sum(a)/sum(a+b)`.
The next command will be perform a sliding windows analysis:
```Bash
/kuacc/users/mbge411/hpc_run/bin/angsd/misc/realSFS fst stats2 results_fst/CEU.pbs.fst.idx -win 50000 -step 10000 -whichFST 1 > results_fst/CEU.pbs.fst.txt
```
We can take a look at the output file:
```Bash
less -S results_fst/CEU.pbs.txt
```
The header is:
`region  chr     midPos  Nsites  Fst01   Fst02   Fst12   PBS0    PBS1    PBS2`
We are interested in the column `PBS2` which gives the PBS values assuming our population (coded here as 2) being the target population (CEU in this case). Note that negative PBS and FST values are equivalent to 0. We are also provided with the individual FST values. You can see that high values of PBS2 are indeed associated with high values of both Fst02 and Fst12 but not Fst01.
Now let us plot our results along with the gene annotation. We can do this using the following command:
```Bash
Rscript scripts/plotPBS.R results_fst/CEU.pbs.fst.txt results_fst/CEU.pbs.fst.txt.pdf 136 137
```
