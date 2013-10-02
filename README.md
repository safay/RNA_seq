RNA_seq
=======

Scripts related to RNA sequencing, transcriptomics, and analysis of differential gene expression

plotSeqLengthDist.R
	This command-line script makes a plot of sequence length distributions, comparing multiple files, initially built to compare the output from different de novo assemblies of transcriptomes.  See the script header for usage.

xprs_2_edgeR.R
	This script is designed to do a series of analyses to explore differential gene expression from an RNA-seq dataset, particularly RNA-seq data from non-model organisms with only a de-novo assembly.  It needs two kinds of files as input: 
	1) a set of xprs.results files from eXpress, one from each library (http://bio.math.berkeley.edu/eXpress/overview.html) 
	2) an annotation file from Trinotate (http://trinotate.sourceforge.net/).  
The script currently is meant to be used interactively in an R IDE like R Studio.
