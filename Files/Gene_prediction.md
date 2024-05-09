# Gene prediction example from human data 


### Objectives:
- Determine location of genes in the sequence with accession number U30787
- Predict genes using [Genscan](http://hollywood.mit.edu/GENSCAN.html)
- Predict genes using [Augustus](http://bioinf.uni-greifswald.de/augustus/submission)

### Find the gene structure information in a genome databases:
- Get sequence from NCBI, in <ins>fasta</ins>  and <ins>genbank</ins> formats (genbank format contains the correct gene structure).
- Locate the sequence on the genome by going to [Ensembl](https://www.ensembl.org/index.html), selecting human, and choosing the blast link. Paste the sequence in the box. You can use BLAT as the search tool.
- The sequence has several differences between the U30787 sequence and the genome reference sequence, which is reflected in a "broken up" set of hits. Select the largest hit in "Alignment summary", note its position in the query (150-4510) and in the genome (chr1:45011305-45015669). Check the alignment to see that there are indels in this large hit.
- Click the corresponding Contigview to view the genome sequence in the matching region with annotation. If you display gene information, you will see that the URO-D gene maps on the sequence, and what the gene structure looks like. You will find many isoforms coming from the same region.
- In the contig view, the blast hits are shown. By default it shows the selected hit with 2kb flanking regions. To see the hit exactly, remove these flanking regions (change location to 1:45010950-45015575, and click "Go"). 
- Use "Export Data" in the left menu and select Genbank or EMBL format to download the data in flat file format. You can easily locate features in this file, but remember that these locations are according to the genome reference sequence! As our query sequence differs (hit starts at 150, indels are present in the hit), you cannot simply transfer these data to our sequence.
- You can find the matching isoform (to the one given in the U30787 genbank file) by searching for the protein_id AAC50482 given in the genbank file in the generated genbank file: it is a db_xref of the isoform identified by transcript_id ENST00000246337 or RefSeq_mRNA NM_000374.

### Predict genes using genscan
- Connect to the [genscan](http://hollywood.mit.edu/GENSCAN.html) server.
- Select organism (vertebrate).
- Paste the DNA sequence in RAW format. Do not copy the fasta header line, only the sequence.
Run GENSCAN.
- Compare the prediction with the annotation in the genbank file. The results contain mostly features of type "Intr", which may be confusing: This does NOT stand for intron, but means that an "internal exon" has been found. Almost all exons are correct. Only the first mostly non-coding exon 1107...1126 is missing, with a wrong prediction at 739..851 instead, and the terminal exon ends where the CDS ends. (If these do not match, you have probably forgotten to remove the fasta header line)

### Predict genes using Augustus
- Connect to the [Augustus](http://bioinf.uni-greifswald.de/augustus/submission) server.
- Select organism (Homo sapiens). Keep "both strands" and "predict any number of (possibly partial) genes". You could set "Alternative transcripts" to "many" if you want to see many potential isoforms. Keep the default here (few), as we do not know any isoforms from the genbank file anyway.
- Paste the DNA sequence in fasta format. WARNING, if you copy the example fasta file as is, you will get erroneous results. Give it a simple fasta name, e.g. ">test"
- Run AUGUSTUS.
- Augustus returns results in gff format. Compare the prediction with the real annotation in the genbank file. The prediction of the second isoform is mostly the same (except start and end, and one missing exon:3279..3416).

