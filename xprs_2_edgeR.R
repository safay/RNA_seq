#!/usr/bin/R

# R script using edgeR for analyzing differential gene expression using count data from eXpress and annotations from Trinotate
# by Scott Fay 
# last updated 1 Oct 2013
# This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 Unported License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/3.0/.

# install needed packages:
#install.packages("ggplot2")
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("edgeR")

library(ggplot2)
library(edgeR)


##########
#
# Load in files with expression estimates (e.g., from eXpress) and FPKM values
# This part needs to be customized for your particular analysis
# ensure you use tot_counts for raw count data in edgeR
#
##########

# set working directory
# this should be a directory with a folder with each library's eXpress output in it, "results.xprs"
setwd("/Volumes/KINGSTON/RKC_DGE_xprs_only_wkshp")
# set directory where you want this script to save its output
out_dir <- "/Volumes/KINGSTON/DGE_out/"

# load in eXpress output files
A4 <- read.delim('./RKC2_A4_index1_ATCACG_L002_xprs_out/results.xprs')
A5 <- read.delim('./RKC2_A5_index25_ACTGAT_L002_xprs_out/results.xprs')
A6 <- read.delim('./RKC2_A6_index27_ATTCCT_L002_xprs_out/results.xprs')
A7 <- read.delim('./RKC2_A7_index1_ATCACG_L005_xprs_out/results.xprs')
A8 <- read.delim('./RKC2_A8_index13_AGTCAA_L002_xprs_out/results.xprs')
B4 <- read.delim('./RKC2_B4_index2_CGATGT_L002_xprs_out/results.xprs')
B5 <- read.delim('./RKC2_B5_index6_GCCAAT_L002_xprs_out/results.xprs')
B6 <- read.delim('./RKC2_B6_index14_AGTTCC_L002_xprs_out/results.xprs')
B7 <- read.delim('./RKC2_B7_index8_ACTTGA_L005_xprs_out/results.xprs')
B8 <- read.delim('./RKC2_B8_index16_CCGTCC_L005_xprs_out/results.xprs')
C1 <- read.delim('./RKC2_C1_index15_ATGTCA_L002_xprs_out/results.xprs')
C2 <- read.delim('./RKC2_C2_index7_CAGATC_L002_xprs_out/results.xprs')
C3 <- read.delim('./RKC2_C3_index20_GTGGCC_L002_xprs_out/results.xprs')
C4 <- read.delim('./RKC2_C4_index9_GATCAG_L005_xprs_out/results.xprs')
C5 <- read.delim('./RKC2_C5_index18_GTCCGC_L005_xprs_out/results.xprs')
D1 <- read.delim('./RKC2_D1_index16_CCGTCC_L002_xprs_out/results.xprs')
D2 <- read.delim('./RKC2_D2_index8_ACTTGA_L002_xprs_out/results.xprs')
D3 <- read.delim('./RKC2_D3_index21_GTTTCG_L002_xprs_out/results.xprs')
D4 <- read.delim('./RKC2_D4_index10_TAGCTT_L005_xprs_out/results.xprs')
D5 <- read.delim('./RKC2_D5_index19_GTGAAA_L005_xprs_out/results.xprs')
E1 <- read.delim('./RKC2_E1_index19_GTGAAA_L002_xprs_out/results.xprs')
E2 <- read.delim('./RKC2_E2_index9_GATCAG_L002_xprs_out/results.xprs')
E3 <- read.delim('./RKC2_E3_index3_TTAGGC_L005_xprs_out/results.xprs')
E4 <- read.delim('./RKC2_E4_index11_GGCTAC_L005_xprs_out/results.xprs')
E5 <- read.delim('./RKC2_E5_index20_GTGGCC_L005_xprs_out/results.xprs')
F1 <- read.delim('./RKC2_F1_index18_GTCCGC_L002_xprs_out/results.xprs')
F2 <- read.delim('./RKC2_F2_index10_TAGCTT_L002_xprs_out/results.xprs')
F3 <- read.delim('./RKC2_F3_index4_TGACCA_L005_xprs_out/results.xprs')
F4 <- read.delim('./RKC2_F4_index12_CTTGTA_L005_xprs_out/results.xprs')
F5 <- read.delim('./RKC2_F5_index21_GTTTCG_L005_xprs_out/results.xprs')
G1 <- read.delim('./RKC2_G1_index3_TTAGGC_L002_xprs_out/results.xprs')
G2 <- read.delim('./RKC2_G2_index11_GGCTAC_L002_xprs_out/results.xprs')
G3 <- read.delim('./RKC2_G3_index5_ACAGTG_L005_xprs_out/results.xprs')
G4 <- read.delim('./RKC2_G4_index13_AGTCAA_L005_xprs_out/results.xprs')
G5 <- read.delim('./RKC2_G5_index22_CGTACG_L005_xprs_out/results.xprs')
H1 <- read.delim('./RKC2_H1_index4_TGACCA_L002_xprs_out/results.xprs')
H2 <- read.delim('./RKC2_H2_index12_CTTGTA_L002_xprs_out/results.xprs')
H3 <- read.delim('./RKC2_H3_index6_GCCAAT_L005_xprs_out/results.xprs')
H4 <- read.delim('./RKC2_H4_index14_AGTTCC_L005_xprs_out/results.xprs')
H5 <- read.delim('./RKC2_H5_index23_GAGTGG_L005_xprs_out/results.xprs')
J1 <- read.delim('./RKC2_J1_index5_ACAGTG_L002_xprs_out/results.xprs')
J2 <- read.delim('./RKC2_J2_index2_CGATGT_L005_xprs_out/results.xprs')
J3 <- read.delim('./RKC2_J3_index7_CAGATC_L005_xprs_out/results.xprs')
J4 <- read.delim('./RKC2_J4_index15_ATGTCA_L005_xprs_out/results.xprs')
J5 <- read.delim('./RKC2_J5_index25_ACTGAT_L005_xprs_out/results.xprs')

# make fpkm data frames
A4_fpkm <- data.frame(target_id=A4$target_id, A4_fpkm=A4$fpkm)
A5_fpkm <- data.frame(target_id=A5$target_id, A5_fpkm=A5$fpkm)
A6_fpkm <- data.frame(target_id=A6$target_id, A6_fpkm=A6$fpkm)
A7_fpkm <- data.frame(target_id=A7$target_id, A7_fpkm=A7$fpkm)
A8_fpkm <- data.frame(target_id=A8$target_id, A8_fpkm=A8$fpkm)
B4_fpkm <- data.frame(target_id=B4$target_id, B4_fpkm=B4$fpkm)
B5_fpkm <- data.frame(target_id=B5$target_id, B5_fpkm=B5$fpkm)
B6_fpkm <- data.frame(target_id=B6$target_id, B6_fpkm=B6$fpkm)
B7_fpkm <- data.frame(target_id=B7$target_id, B7_fpkm=B7$fpkm)
B8_fpkm <- data.frame(target_id=B8$target_id, B8_fpkm=B8$fpkm)
C1_fpkm <- data.frame(target_id=C1$target_id, C1_fpkm=C1$fpkm)
C2_fpkm <- data.frame(target_id=C2$target_id, C2_fpkm=C2$fpkm)
C3_fpkm <- data.frame(target_id=C3$target_id, C3_fpkm=C3$fpkm)
C4_fpkm <- data.frame(target_id=C4$target_id, C4_fpkm=C4$fpkm)
C5_fpkm <- data.frame(target_id=C5$target_id, C5_fpkm=C5$fpkm)
D1_fpkm <- data.frame(target_id=D1$target_id, D1_fpkm=D1$fpkm)
D2_fpkm <- data.frame(target_id=D2$target_id, D2_fpkm=D2$fpkm)
D3_fpkm <- data.frame(target_id=D3$target_id, D3_fpkm=D3$fpkm)
D4_fpkm <- data.frame(target_id=D4$target_id, D4_fpkm=D4$fpkm)
D5_fpkm <- data.frame(target_id=D5$target_id, D5_fpkm=D5$fpkm)
E1_fpkm <- data.frame(target_id=E1$target_id, E1_fpkm=E1$fpkm)
E2_fpkm <- data.frame(target_id=E2$target_id, E2_fpkm=E2$fpkm)
E3_fpkm <- data.frame(target_id=E3$target_id, E3_fpkm=E3$fpkm)
E4_fpkm <- data.frame(target_id=E4$target_id, E4_fpkm=E4$fpkm)
E5_fpkm <- data.frame(target_id=E5$target_id, E5_fpkm=E5$fpkm)
F1_fpkm <- data.frame(target_id=F1$target_id, F1_fpkm=F1$fpkm)
F2_fpkm <- data.frame(target_id=F2$target_id, F2_fpkm=F2$fpkm)
F3_fpkm <- data.frame(target_id=F3$target_id, F3_fpkm=F3$fpkm)
F4_fpkm <- data.frame(target_id=F4$target_id, F4_fpkm=F4$fpkm)
F5_fpkm <- data.frame(target_id=F5$target_id, F5_fpkm=F5$fpkm)
G1_fpkm <- data.frame(target_id=G1$target_id, G1_fpkm=G1$fpkm)
G2_fpkm <- data.frame(target_id=G2$target_id, G2_fpkm=G2$fpkm)
G3_fpkm <- data.frame(target_id=G3$target_id, G3_fpkm=G3$fpkm)
G4_fpkm <- data.frame(target_id=G4$target_id, G4_fpkm=G4$fpkm)
G5_fpkm <- data.frame(target_id=G5$target_id, G5_fpkm=G5$fpkm)
H1_fpkm <- data.frame(target_id=H1$target_id, H1_fpkm=H1$fpkm)
H2_fpkm <- data.frame(target_id=H2$target_id, H2_fpkm=H2$fpkm)
H3_fpkm <- data.frame(target_id=H3$target_id, H3_fpkm=H3$fpkm)
H4_fpkm <- data.frame(target_id=H4$target_id, H4_fpkm=H4$fpkm)
H5_fpkm <- data.frame(target_id=H5$target_id, H5_fpkm=H5$fpkm)
J1_fpkm <- data.frame(target_id=J1$target_id, J1_fpkm=J1$fpkm)
J2_fpkm <- data.frame(target_id=J2$target_id, J2_fpkm=J2$fpkm)
J3_fpkm <- data.frame(target_id=J3$target_id, J3_fpkm=J3$fpkm)
J4_fpkm <- data.frame(target_id=J4$target_id, J4_fpkm=J4$fpkm)
J5_fpkm <- data.frame(target_id=J5$target_id, J5_fpkm=J5$fpkm)

# merge fpkm data frames into one new data frame
all_fpkm <- merge(A4_fpkm, A5_fpkm)
all_fpkm <- merge(all_fpkm, A6_fpkm)
all_fpkm <- merge(all_fpkm, A7_fpkm)
all_fpkm <- merge(all_fpkm, A8_fpkm)
all_fpkm <- merge(all_fpkm, B4_fpkm)
all_fpkm <- merge(all_fpkm, B5_fpkm)
all_fpkm <- merge(all_fpkm, B6_fpkm)
all_fpkm <- merge(all_fpkm, B7_fpkm)
all_fpkm <- merge(all_fpkm, B8_fpkm)
all_fpkm <- merge(all_fpkm, C1_fpkm)
all_fpkm <- merge(all_fpkm, C2_fpkm)
all_fpkm <- merge(all_fpkm, C3_fpkm)
all_fpkm <- merge(all_fpkm, C4_fpkm)
all_fpkm <- merge(all_fpkm, C5_fpkm)
all_fpkm <- merge(all_fpkm, D1_fpkm)
all_fpkm <- merge(all_fpkm, D2_fpkm)
all_fpkm <- merge(all_fpkm, D3_fpkm)
all_fpkm <- merge(all_fpkm, D4_fpkm)
all_fpkm <- merge(all_fpkm, D5_fpkm)
all_fpkm <- merge(all_fpkm, E1_fpkm)
all_fpkm <- merge(all_fpkm, E2_fpkm)
all_fpkm <- merge(all_fpkm, E3_fpkm)
all_fpkm <- merge(all_fpkm, E4_fpkm)
all_fpkm <- merge(all_fpkm, E5_fpkm)
all_fpkm <- merge(all_fpkm, F1_fpkm)
all_fpkm <- merge(all_fpkm, F2_fpkm)
all_fpkm <- merge(all_fpkm, F3_fpkm)
all_fpkm <- merge(all_fpkm, F4_fpkm)
all_fpkm <- merge(all_fpkm, F5_fpkm)
all_fpkm <- merge(all_fpkm, G1_fpkm)
all_fpkm <- merge(all_fpkm, G2_fpkm)
all_fpkm <- merge(all_fpkm, G3_fpkm)
all_fpkm <- merge(all_fpkm, G4_fpkm)
all_fpkm <- merge(all_fpkm, G5_fpkm)
all_fpkm <- merge(all_fpkm, H1_fpkm)
all_fpkm <- merge(all_fpkm, H2_fpkm)
all_fpkm <- merge(all_fpkm, H3_fpkm)
all_fpkm <- merge(all_fpkm, H4_fpkm)
all_fpkm <- merge(all_fpkm, H5_fpkm)
all_fpkm <- merge(all_fpkm, J1_fpkm)
all_fpkm <- merge(all_fpkm, J2_fpkm)
all_fpkm <- merge(all_fpkm, J3_fpkm)
all_fpkm <- merge(all_fpkm, J4_fpkm)
all_fpkm <- merge(all_fpkm, J5_fpkm)

# make total count data frames
A4_count <- data.frame(target_id=A4$target_id, A4_count=A4$tot_counts)
A5_count <- data.frame(target_id=A5$target_id, A5_count=A5$tot_counts)
A6_count <- data.frame(target_id=A6$target_id, A6_count=A6$tot_counts)
A7_count <- data.frame(target_id=A7$target_id, A7_count=A7$tot_counts)
A8_count <- data.frame(target_id=A8$target_id, A8_count=A8$tot_counts)
B4_count <- data.frame(target_id=B4$target_id, B4_count=B4$tot_counts)
B5_count <- data.frame(target_id=B5$target_id, B5_count=B5$tot_counts)
B6_count <- data.frame(target_id=B6$target_id, B6_count=B6$tot_counts)
B7_count <- data.frame(target_id=B7$target_id, B7_count=B7$tot_counts)
B8_count <- data.frame(target_id=B8$target_id, B8_count=B8$tot_counts)
C1_count <- data.frame(target_id=C1$target_id, C1_count=C1$tot_counts)
C2_count <- data.frame(target_id=C2$target_id, C2_count=C2$tot_counts)
C3_count <- data.frame(target_id=C3$target_id, C3_count=C3$tot_counts)
C4_count <- data.frame(target_id=C4$target_id, C4_count=C4$tot_counts)
C5_count <- data.frame(target_id=C5$target_id, C5_count=C5$tot_counts)
D1_count <- data.frame(target_id=D1$target_id, D1_count=D1$tot_counts)
D2_count <- data.frame(target_id=D2$target_id, D2_count=D2$tot_counts)
D3_count <- data.frame(target_id=D3$target_id, D3_count=D3$tot_counts)
D4_count <- data.frame(target_id=D4$target_id, D4_count=D4$tot_counts)
D5_count <- data.frame(target_id=D5$target_id, D5_count=D5$tot_counts)
E1_count <- data.frame(target_id=E1$target_id, E1_count=E1$tot_counts)
E2_count <- data.frame(target_id=E2$target_id, E2_count=E2$tot_counts)
E3_count <- data.frame(target_id=E3$target_id, E3_count=E3$tot_counts)
E4_count <- data.frame(target_id=E4$target_id, E4_count=E4$tot_counts)
E5_count <- data.frame(target_id=E5$target_id, E5_count=E5$tot_counts)
F1_count <- data.frame(target_id=F1$target_id, F1_count=F1$tot_counts)
F2_count <- data.frame(target_id=F2$target_id, F2_count=F2$tot_counts)
F3_count <- data.frame(target_id=F3$target_id, F3_count=F3$tot_counts)
F4_count <- data.frame(target_id=F4$target_id, F4_count=F4$tot_counts)
F5_count <- data.frame(target_id=F5$target_id, F5_count=F5$tot_counts)
G1_count <- data.frame(target_id=G1$target_id, G1_count=G1$tot_counts)
G2_count <- data.frame(target_id=G2$target_id, G2_count=G2$tot_counts)
G3_count <- data.frame(target_id=G3$target_id, G3_count=G3$tot_counts)
G4_count <- data.frame(target_id=G4$target_id, G4_count=G4$tot_counts)
G5_count <- data.frame(target_id=G5$target_id, G5_count=G5$tot_counts)
H1_count <- data.frame(target_id=H1$target_id, H1_count=H1$tot_counts)
H2_count <- data.frame(target_id=H2$target_id, H2_count=H2$tot_counts)
H3_count <- data.frame(target_id=H3$target_id, H3_count=H3$tot_counts)
H4_count <- data.frame(target_id=H4$target_id, H4_count=H4$tot_counts)
H5_count <- data.frame(target_id=H5$target_id, H5_count=H5$tot_counts)
J1_count <- data.frame(target_id=J1$target_id, J1_count=J1$tot_counts)
J2_count <- data.frame(target_id=J2$target_id, J2_count=J2$tot_counts)
J3_count <- data.frame(target_id=J3$target_id, J3_count=J3$tot_counts)
J4_count <- data.frame(target_id=J4$target_id, J4_count=J4$tot_counts)
J5_count <- data.frame(target_id=J5$target_id, J5_count=J5$tot_counts)

# merge total counts into one data frame
all_count <- merge(A4_count, A5_count)
all_count <- merge(all_count, A6_count)
all_count <- merge(all_count, A7_count)
all_count <- merge(all_count, A8_count)
all_count <- merge(all_count, B4_count)
all_count <- merge(all_count, B5_count)
all_count <- merge(all_count, B6_count)
all_count <- merge(all_count, B7_count)
all_count <- merge(all_count, B8_count)
all_count <- merge(all_count, C1_count)
all_count <- merge(all_count, C2_count)
all_count <- merge(all_count, C3_count)
all_count <- merge(all_count, C4_count)
all_count <- merge(all_count, C5_count)
all_count <- merge(all_count, D1_count)
all_count <- merge(all_count, D2_count)
all_count <- merge(all_count, D3_count)
all_count <- merge(all_count, D4_count)
all_count <- merge(all_count, D5_count)
all_count <- merge(all_count, E1_count)
all_count <- merge(all_count, E2_count)
all_count <- merge(all_count, E3_count)
all_count <- merge(all_count, E4_count)
all_count <- merge(all_count, E5_count)
all_count <- merge(all_count, F1_count)
all_count <- merge(all_count, F2_count)
all_count <- merge(all_count, F3_count)
all_count <- merge(all_count, F4_count)
all_count <- merge(all_count, F5_count)
all_count <- merge(all_count, G1_count)
all_count <- merge(all_count, G2_count)
all_count <- merge(all_count, G3_count)
all_count <- merge(all_count, G4_count)
all_count <- merge(all_count, G5_count)
all_count <- merge(all_count, H1_count)
all_count <- merge(all_count, H2_count)
all_count <- merge(all_count, H3_count)
all_count <- merge(all_count, H4_count)
all_count <- merge(all_count, H5_count)
all_count <- merge(all_count, J1_count)
all_count <- merge(all_count, J2_count)
all_count <- merge(all_count, J3_count)
all_count <- merge(all_count, J4_count)
all_count <- merge(all_count, J5_count)

# create row names from target_id
rownames(all_count) <- all_count$target_id
rownames(all_fpkm) <- all_fpkm$target_id
# remove target_id column
all_count$target_id <- NULL
all_fpkm$target_id <- NULL
# now the data is ready

# define groups
pH <- factor(c(rep("amb",5), rep("low",5), rep("amb",5), rep("low",5), rep("mid",5), rep("amb",5), rep("mid",5), rep("low",5), rep("mid",5)))
temp <- factor(c(rep("13C", 5), rep("11C", 5), rep("11C",5), rep("14C", 5), rep("13C",5), rep("14C", 5), rep("14C", 5), rep("13C", 5), rep("11C",5)))
group <- paste(pH,temp, sep="")
pH
temp
group

# remove spike option
# use this if you have certain contaminating sequences you want to remove
spikes <- scan("~/path/to/spike_seqIDs.txt", what = "character")
# remove the spikes transcripts from the tables
all_count <- all_count[!rownames(all_count) %in% spikes, ]
all_fpkm <- all_fpkm[!rownames(all_fpkm) %in% spikes, ]

# Get annotations
# "trinotate_annotation_report.txt" is a Trinotate output excel file exported in excel as tab-delimited txt
transcripts <- read.delim("/path/to/trinotate_annotation_report.txt")

##########
#
# Set up differential gene expression (DGE) analysis, generate summary plots:
#   define a model matrix
#   estimate dispersion
#   make a plot of biological coefficient of variation vs. log(CPM)
#   make a multidimensional scaling plot of total counts
#   fit GLM for each transcript
#
##########


# DGEList makes an EdgeR object
y <- DGEList(counts=all_count, group=group)
# define the statistical model, a design matrix using the model.matrix function
# crossed design, same as (~ pH + temp + pH:temp ... see section 3.3.4 in edgeR documentation)
design <- model.matrix(~pH*temp, data=y$samples) 
# view the design matrix
design
colnames(design)
# filter out transcripts that do not have at least two reads out of 1,000,000 reads mapped in at least 4 samples
y <- y[rowSums(1e+06 * y$counts/expandAsMatrix(y$samples$lib.size, dim(y)) > 2) >= 4, ]
# reset the library sizes after filtering
y$samples$lib.size <- colSums(y$counts)
# TMM normalization compensates not just for library size but also the relative expression level among transcripts
y <- calcNormFactors(y)
# estimate dispersion for GLM fit, common, trended, and tagwise
# see edgeR docs for references on dispersion estimates, e.g., "?estimateGLMCommonDisp"
y <- estimateGLMCommonDisp(y,design) 
# estimate trended dispersion for use in tagwise dispersion estimate
y <- estimateGLMTrendedDisp(y,design)
# estimate tagwise dispersion to be used in glmFit()
y <- estimateGLMTagwiseDisp(y,design)
# plot genewise biological coefficient of variation against gene abundance
plotBCV(y)
# make plot as a pdf
pdf(file=paste(out_dir, "BCV_plot.pdf", sep=""), height=6, width=6)
plotBCV(y, main = "Biological Coefficient of Variation")
dev.off()
# MDS plot
plotMDS(y , main = "MDS Plot for Count Data", labels = colnames( y$counts ), cex=0.7)
# make a pdf
pdf(file=paste(out_dir, "MDS_plot.pdf", sep=""), height=6, width=6)
plotMDS(y , main = "MDS Plot for Count Data", labels = colnames( y$counts ), cex=0.7)
dev.off()

# GLM
# fit the negative binomial generalized linear model (GLM) for each tag, creating a new fit object
fit <- glmFit(y, design)
colnames(fit)

#####
#
# Find significantly differentially expressed genes between comparisons
#
#####

#####
# Function: get_DEGs
# perform LRT (likelihood ratio test) on certain factors
# returns a data frame that contains both toptags and annotation data
# arguments:
#   compare = coefficient for comparison
#   annot = transcripts annotation data frame
#   fit = glmFit object (fit, in this case)
#   fdr = false discovery rate (default 0.05)
#####

get_DEGs <- function(compare, annot, fit, fdr=0.05) {
  lrt <- glmLRT(fit, coef=compare)
  DEG <- summary(decideTestsDGE(lrt, p=fdr, adjust="BH")) # gives numbers of genes at FDR < pval
  Num_DEG <- (DEG[1] + DEG[-1])[2] # number of DEG based on fdr
  tTags <- topTags(lrt, n=Num_DEG) # get list of DEGs, "topTags"
  cat("Comparison:", tTags$comparison, "\n")
  cat("Total number of genes:", nrow(lrt$table), "\n")
  cat("Number of differentially expressed genes with FDR <", fdr, "=", nrow(tTags), "\n")
  DEG
  cat("...saving an MAplot:", paste(out_dir,"MAplot_", tTags$comparison, "_FDR_", fdr, ".pdf", sep=""), "\n" )
  pdf(file=paste(out_dir,"MAplot_", tTags$comparison, "_FDR_", fdr, ".pdf", sep=""), height=6, width=6)
  detags <- rownames(tTags$table) # n= has to be the number of significantly differentially expressed genes
  plotSmear(lrt,de.tags=detags, cex=0.5) # plot of fold change given CPM, red for those < FDR
  abline(h=c(-1,1), col="dodgerblue") # blue line at twofold change
  dev.off()
  # Merge annotations with DGE stats
  tTags_frame <- tTags$table # make data frame from tTags
  tTags_frame$transcript_id <- rownames(tTags_frame)
  join <- merge(tTags_frame, annot) # gets the intersection of detags and transcripts, i.e., the annotations for the DE transcripts
  cat("percent all transcripts with ORFs:", sum(annot$prot_id != ".") * 100 / nrow(annot), "\n")
  cat("percent toptags with ORFs:", sum(join$prot_id != ".") * 100 / nrow(join), "\n")
  cat("percent all transcripts with Pfam or TopBlastHit annotation:", sum((annot$Pfam != ".") | (annot$TopBlastHit != ".")) * 100 / nrow(annot), "\n")
  cat("percent toptags with Pfam or TopBlastHit annotation:", sum((join$Pfam != ".") | (join$TopBlastHit != ".")) * 100 / nrow(join), "\n")
  cat("number of genes upregulated, (FC > 2):", sum(join$logFC > 1), "\n")
  cat("number of genes downregulated, (FC < 0.5):", sum(join$logFC < -1 ), "\n")
  cat("...saving .csv file with all topTags:", paste( out_dir, "topTags_", tTags$comparison, "_", fdr, "_FDR.csv", sep="" ), "\n" )
  write.csv( join, file=paste(out_dir, "topTags_", tTags$comparison, "_", fdr, "_FDR.csv", sep="" ) )
  cat("...saving .csv file with pFam OR TopBlastHit annotated topTags:", paste(out_dir, "topTags_", tTags$comparison, "_", fdr, "_FDR_only_w_annot.csv", sep=""), "\n" )
  write.csv( join[ (join$Pfam != ".") | (join$TopBlastHit != ".") , ], file=paste(out_dir, "topTags_", tTags$comparison, "_", fdr, "_FDR_only_w_PFamAnnot.csv", sep="" ) )
  join
}

# call the get_DEGs() function for each comparison
DEGlist_pHlow <- get_DEGs(compare=2, annot=transcripts, fit=fit, fdr=0.00001)
DEGlist_pHmid <- get_DEGs(compare=3, annot=transcripts, fit=fit, fdr=0.00001)
DEGlist_temp13C <- get_DEGs(compare=4, annot=transcripts, fit=fit, fdr=0.00001)
DEGlist_temp14C <- get_DEGs(compare=5, annot=transcripts, fit=fit, fdr=0.00001)
DEGlist_pHlow_temp13C <- get_DEGs(compare=6, annot=transcripts, fit=fit, fdr=0.00001)
DEGlist_pHmid_temp13C <- get_DEGs(compare=7, annot=transcripts, fit=fit, fdr=0.00001)
DEGlist_pHlow_temp14C <- get_DEGs(compare=8, annot=transcripts, fit=fit, fdr=0.00001)
DEGlist_pHmid_temp14C <- get_DEGs(compare=9, annot=transcripts, fit=fit, fdr=0.00001)

# sort by logFC
sorted_tags <- annotated_DEGs_pHlow[with(annotated_DEGs_pHlow, order(logFC)), ]

# write a .csv based on a threshold FC, in this case, ( FC < 0.5 ) and ( FC > 2 )
write.csv(annotated_DEGs_pHlow[ abs(annotated_DEGs_pHlow$logFC) > 1 , ], file=paste( out_dir, "topTags_DEGlist_pHlow_FDR0_00001_PfamHits_logFC_greater_than_2.csv", sep="" )

#####
# Function: deg_heatmap
#   fpkm_frame, the dataframe made above, "all_fpkm"
#   tTagList, a list yeilded from the get_DGEs function
#   crit_mean_fpkm, critical value for filtering out rows based on mean FPKM across all samples.  Default = 0, i.e., all samples.
#####

DEG_heatmap <- function(fpkm_frame, tTagList, crit_mean_fpkm = 0) {
  filtered_fpkm_frame <- fpkm_frame[ ( apply(fpkm_frame, 1, mean) > crit_mean_fpkm ) , ]  # filter fpkm_frame by rows above a critical mean fpkm value, crit_mean_fpkm
  # reorder by pH:
  filtered_fpkm_frame <- filtered_fpkm_frame[c("C1_fpkm", "C2_fpkm", "C3_fpkm", "C4_fpkm", "C5_fpkm", "A4_fpkm", "A5_fpkm", "A6_fpkm", "A7_fpkm", "A8_fpkm", "F1_fpkm", "F2_fpkm", "F3_fpkm", "F4_fpkm", "F5_fpkm", "J1_fpkm", "J2_fpkm", "J3_fpkm", "J4_fpkm", "J5_fpkm", "E1_fpkm", "E2_fpkm", "E3_fpkm", "E4_fpkm", "E5_fpkm", "G1_fpkm", "G2_fpkm", "G3_fpkm", "G4_fpkm", "G5_fpkm", "B4_fpkm", "B5_fpkm", "B6_fpkm", "B7_fpkm", "B8_fpkm", "H1_fpkm", "H2_fpkm", "H3_fpkm", "H4_fpkm", "H5_fpkm", "D1_fpkm", "D2_fpkm", "D3_fpkm", "D4_fpkm", "D5_fpkm" ) ]
  fpkm_matrix <- data.matrix(filtered_fpkm_frame[ ( rownames(filtered_fpkm_frame) %in% tTagList$transcript_id ) ,  ]) # make matrix of only those transcripts found in tTagList
  n=256
  heatmap(fpkm_matrix, Colv=NA, col = rainbow(n, s = 1, v = 1, start = 0, end = max(1,n - 1)/n, alpha = 1), scale="column", margins=c(5,10), cexRow = 0.5, cexCol = 0.7)
}

# Example usage: 
DEG_heatmap(all_fpkm, DEGlist_no_interaction_pHlow)
DEG_heatmap(all_fpkm, DEGlist_no_interaction_pHlow, 100)

##########
#
# Functions for genewise data analysis
#
##########

####
# Function: get_gene_data
# examine expression of an individual gene based on transcript_id
####

get_gene_data <- function(gene, fpkm) {
  A_mean <- apply(fpkm[gene,1:5], 1, mean)
  B_mean <- apply(fpkm[gene,6:10], 1, mean)
  C_mean <- apply(fpkm[gene,11:15], 1, mean)
  D_mean <- apply(fpkm[gene,16:20], 1, mean)
  E_mean <- apply(fpkm[gene,21:25], 1, mean)
  F_mean <- apply(fpkm[gene,26:30], 1, mean)
  G_mean <- apply(fpkm[gene,31:35], 1, mean)
  H_mean <- apply(fpkm[gene,36:40], 1, mean)
  J_mean <- apply(fpkm[gene,41:45], 1, mean)
  A_sd <- apply(fpkm[gene,1:5], 1, sd)
  B_sd <- apply(fpkm[gene,6:10], 1, sd)
  C_sd <- apply(fpkm[gene,11:15], 1, sd)
  D_sd <- apply(fpkm[gene,16:20], 1, sd)
  E_sd <- apply(fpkm[gene,21:25], 1, sd)
  F_sd <- apply(fpkm[gene,26:30], 1, sd)
  G_sd <- apply(fpkm[gene,31:35], 1, sd)
  H_sd <- apply(fpkm[gene,36:40], 1, sd)
  J_sd <- apply(fpkm[gene,41:45], 1, sd)
  fpkm_mean <- c(A_mean, B_mean, C_mean, D_mean, E_mean, F_mean, G_mean, H_mean, J_mean)
  SD <- c(A_sd, B_sd, C_sd, D_sd, E_sd, F_sd, G_sd, H_sd, J_sd)
  names <- c("A", "B", "C", "D", "E", "F", "G", "H", "J")
  pH_grp <- c("amb", "low", "amb", "low", "mid", "amb", "mid", "low", "mid")
  temp_grp <- c("13C", "11C", "11C", "14C", "13C", "14C", "14C", "13C", "11C")
  out_frame <- data.frame(fpkm_mean, SD, pH_grp, temp_grp, names)
  out_frame$pH_grp <- factor( out_frame$pH_grp, levels = c( "amb", "mid", "low" ) )
  return(out_frame)
}

# get fpkm data for a specific gene
#gene_data <- get_gene_data("comp106087_c0_seq1", all_fpkm)

#####
# Function: get_annot
# gets annotation for a particular gene
#####

get_annot <- function(gene_in, annotation) {
  L <- annotation$transcript_id == gene_in
  cat("transcript: ", as.vector(annotation[L,]$transcript_id), "\n")
  cat("top blast hit: ", as.vector(annotation[L,]$TopBlastHit), "\n")
  cat("gene ontology: ", as.vector(annotation[L,]$gene_ontology), "\n")
  cat("PFam: ", as.vector(annotation[L,]$Pfam), "\n")
  cat("Signal Peptide: ", as.vector(annotation[L,]$SignalP), "\n")
  cat("transmembrane domain: ", as.vector(annotation[L,]$TmHMM), "\n")
  cat("eggnog: ", as.vector(annotation[L,]$eggnog), "\n")
}

# test getting annotation data
#get_annot("comp99839_c0_seq1", transcripts)

#####
# Function: plot_genewise_fpkm
# plots fpkm for the specified gene
#####

plot_genewise_fpkm <- function(gene_in, plot_title, add_title) {
  k <- ggplot(gene_in, aes( x = temp_grp, y = fpkm_mean, ymax = fpkm_mean + SD, ymin = fpkm_mean - SD, fill = temp_grp ))
  k + geom_bar(stat="identity") + facet_grid( . ~ pH_grp ) + geom_errorbar( width = 0.3) + scale_fill_brewer(type="qual", palette="Paired") + labs(title=paste(plot_title, add_title))
  k + geom_bar(stat="identity") + facet_grid( . ~ pH_grp ) + geom_errorbar( width = 0.3) + scale_fill_brewer(type="qual", palette="Paired") + labs(title=paste(plot_title, add_title))
  ggsave(filename = paste(out_dir, plot_title, "_", add_title, "_genewise_FPKM", ".pdf", sep =""), width = 6, height = 6)
  }

# test plotting fpkm data
#plot_genewise_fpkm(gene_data)

#####
# function: plot_interaction
# makes an interaction plot for a given gene
#####

plot_interaction <- function(gene_data, plot_title, add_title) {
  h <- ggplot(gene_data, aes(x = pH_grp, y = fpkm_mean, ymax = fpkm_mean + SD, ymin = fpkm_mean - SD, color=temp_grp))
  h + geom_errorbar(color="black", width=0.1) + geom_point(aes(shape=temp_grp), color="black",size = 6) + geom_line(aes(group=temp_grp), size=0.5, linetype=1) + labs(title=paste(plot_title, add_title), x = "pH treatment", y = "mean FPKM") + scale_color_manual(values=c("#3399FF", "#66CC00", "#CC0066"))
  ggsave(filename = paste(out_dir, plot_title, "_", add_title, "_interactionPlot_FPKM", ".pdf", sep =""), width = 6, height = 6)
  }

#####
# function: GOI_data
# master function to get data and plots from Gene Of Interest
#   gene_name: the "transcript_id" for your gene of interest
#   annotation: the full trinotate annotation report as a data frame
#   fpkm_table: the table of FPKM values 
#   graph_title: a brief gene name you want to use for the graph
#####

GOI_data <- function(gene_name, annotation, fpkm_table, graph_title) {
  gene_data <- get_gene_data(gene_name, fpkm_table)
  get_annot(gene_name, annotation)
  plot_genewise_fpkm(gene_data, gene_name, graph_title)
  plot_interaction(gene_data, gene_name, graph_title)
}

# generate genewise plots based on a gene of interest
GOI_data("comp100107_c0_seq1", transcripts, all_fpkm, "Carbonic Anhydrase")