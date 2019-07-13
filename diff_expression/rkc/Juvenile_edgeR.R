# Scott Fay 
# 17 Sept 2013
# updated June 2014

# first, install packages:
# install.packages("ggplot2")
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("edgeR")

library(ggplot2)
library(edgeR)


##########
#
# Load in files with expression estimates (e.g., from eXpress) and FPKM values
# This part needs to be customized for your particular analysis
#
##########

# set working directory
setwd("/Users/Shared/RKC_project/RKC2/RKC2_June2014/RKC2.1_xprsOut")
# set output directory
out_dir = "/Users/jonathonstillman/Documents/RKC_Transcriptomics/RKC2_Juveniles/"
# below is for working on BerkeleyiMac
#out_dir <- "/Users/Shared/RKC_project/RKC2/RKC2_June2014"

# load in files to get count file for edgeR and FPKM for heatmap etc.
A4 <- read.delim('./RKC2_A4_index1_ATCACG_L002_xprs_out/results.xprs')
A5 <- read.delim('./RKC2_A5_index25_ACTGAT_L002_xprs_out/results.xprs')
A6 <- read.delim('./RKC2_A6_index27_ATTCCT_L002_xprs_out/results.xprs')
C2 <- read.delim('./RKC2_A7_index1_ATCACG_L005_xprs_out/results.xprs')
A8 <- read.delim('./RKC2_A8_index13_AGTCAA_L002_xprs_out/results.xprs')
B4 <- read.delim('./RKC2_B4_index2_CGATGT_L002_xprs_out/results.xprs')
B5 <- read.delim('./RKC2_B5_index6_GCCAAT_L002_xprs_out/results.xprs')
B6 <- read.delim('./RKC2_B6_index14_AGTTCC_L002_xprs_out/results.xprs')
B7 <- read.delim('./RKC2_B7_index8_ACTTGA_L005_xprs_out/results.xprs')
B8 <- read.delim('./RKC2_B8_index16_CCGTCC_L005_xprs_out/results.xprs')
C1 <- read.delim('./RKC2_C1_index15_ATGTCA_L002_xprs_out/results.xprs')
A7 <- read.delim('./RKC2_C2_index7_CAGATC_L002_xprs_out/results.xprs')
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
H1 <- read.delim('./RKC2_G1_index3_TTAGGC_L002_xprs_out/results.xprs')
G2 <- read.delim('./RKC2_G2_index11_GGCTAC_L002_xprs_out/results.xprs')
G3 <- read.delim('./RKC2_G3_index5_ACAGTG_L005_xprs_out/results.xprs')
G4 <- read.delim('./RKC2_G4_index13_AGTCAA_L005_xprs_out/results.xprs')
G5 <- read.delim('./RKC2_G5_index22_CGTACG_L005_xprs_out/results.xprs')
J1 <- read.delim('./RKC2_H1_index4_TGACCA_L002_xprs_out/results.xprs')
H2 <- read.delim('./RKC2_H2_index12_CTTGTA_L002_xprs_out/results.xprs')
H3 <- read.delim('./RKC2_H3_index6_GCCAAT_L005_xprs_out/results.xprs')
H4 <- read.delim('./RKC2_H4_index14_AGTTCC_L005_xprs_out/results.xprs')
H5 <- read.delim('./RKC2_H5_index23_GAGTGG_L005_xprs_out/results.xprs')
G1 <- read.delim('./RKC2_J1_index5_ACAGTG_L002_xprs_out/results.xprs')
J2 <- read.delim('./RKC2_J2_index2_CGATGT_L005_xprs_out/results.xprs')
J3 <- read.delim('./RKC2_J3_index7_CAGATC_L005_xprs_out/results.xprs')
J4 <- read.delim('./RKC2_J4_index15_ATGTCA_L005_xprs_out/results.xprs')
J5 <- read.delim('./RKC2_J5_index25_ACTGAT_L005_xprs_out/results.xprs')

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
pH2 <- factor(c(rep("8.1",5), rep("7.5",5), rep("8.1",5), rep("7.5",5), rep("7.8",5), rep("8.1",5), rep("7.8",5), rep("7.5",5), rep("7.8",5)))
temp <- factor(c(rep("13C", 5), rep("11C", 5), rep("11C",5), rep("14C", 5), rep("13C",5), rep("14C", 5), rep("14C", 5), rep("13C", 5), rep("11C",5)))
group <- paste(pH,temp, sep="")
pH
temp
group

# remove spike option
# retreive the spike IDs generated using blast
spikes <- scan("/Users/scottfay/Google Drive_SF Old/RNA_seq/RKC/RKC_Assemblies", what = "character")
# remove the spikes transcripts from the tables
all_count <- all_count[!rownames(all_count) %in% spikes, ]
all_fpkm <- all_fpkm[!rownames(all_fpkm) %in% spikes, ]

# Get annotations
# "trinotate_annotation_report.txt" is a Trinotate output excel file exported tab-delimited
transcripts <- read.csv("/Users/scottfay/Dropbox/RKC_ALL_Mapped_Annot/ALL_RKC_trinotate_annotation_report.csv", stringsAsFactors=FALSE)

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
design <- model.matrix(~pH*temp, data=y$samples) # crossed design, same as (~ pH + temp + pH:temp ... see section 3.3.4 in edgeR documentation, and double check this meaning with a statistician
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
plotMDS(y , main = "MDS Plot for Count Data", labels = group, cex=0.7)
dev.off()

boxplot(log2(y$counts+1), las=2, main="After filter") # before filter

#####
#
# GLM
#
#####

# fit the negative binomial generalized linear model (GLM) for each tag, creating a new fit object
fit <- glmFit(y, design)
colnames(fit)

# perform likelihood ratio tests for various comparisons, generating LRT object for each

## pH comparisons
lrt_pH_1 <- glmLRT(fit, coef=2)
lrt_pH_2 <- glmLRT(fit, coef=3)
lrt_pH_3 <- glmLRT(fit, contrast=c(0,-1,1,0,0,0,0,0,0))

### temp comparisons
lrt_temp_1 <- glmLRT(fit, coef=4)
lrt_temp_2 <- glmLRT(fit, coef=5)
lrt_temp_3 <- glmLRT(fit, contrast=c(0,0,0,-1,1,0,0,0,0))

### interaction term comparisons
lrt_inter_0 <- glmLRT(fit, coef=6)
lrt_inter_1 <- glmLRT(fit, coef=7)
lrt_inter_2 <- glmLRT(fit, coef=8)
lrt_inter_3 <- glmLRT(fit, coef=9)
lrt_inter_4 <- glmLRT(fit, contrast=c(0,0,0,0,0,-1,1,0,0))
lrt_inter_5 <- glmLRT(fit, contrast=c(0,0,0,0,0,-1,0,1,0))
lrt_inter_6 <- glmLRT(fit, contrast=c(0,0,0,0,0,-1,0,0,1))
lrt_inter_7 <- glmLRT(fit, contrast=c(0,0,0,0,0,0,-1,1,0))
lrt_inter_8 <- glmLRT(fit, contrast=c(0,0,0,0,0,0,-1,0,1))
lrt_inter_9 <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,-1,1))



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
#   lrt = glmLRT object of the comparison you're making
#   annot = transcripts annotation data frame
#   fdr = false discovery rate (default 0.05)
#   critFC = a threshold fold change (default 2-fold)
#   onlyAnnot = flag whether to save only annotated transcripts
#
# saves files:
#   .pdf of the MA plot
#   .csv file of topTags, filtered by critical FC and onlyAnnot flag
#   
#####

get_DEGs <- function(lrt, annot, fdr=0.05, critFC=2, onlyAnnot=FALSE) {
  DEG <- summary(decideTestsDGE(lrt, p=fdr, adjust="BH")) # gives numbers of genes up and downregulated at FDR < pval
  Num_DEG <- (DEG[1] + DEG[-1])[2] # gets total number of DEGs based on the FDR
  tTags <- topTags(lrt, n=Num_DEG) # gets a list of DEGs, "topTags"
  cp <- unlist(strsplit(tTags$comparison, "[* ]"))
  comp <- paste(cp[2], "v", cp[4], sep="") # this and the preceding line make a cleaner text format for the "comparison," used in filenames
  cat("Comparison:", comp, "\n")
  cat("Total number of genes:", nrow(lrt$table), "\n")
  cat("Number of differentially expressed genes with FDR <", fdr, "=", nrow(tTags), "\n")
  # write an MA Plot .pdf
  cat("...saving an MAplot:", paste(out_dir,"MAplot_", comp, "_FDR_", fdr, ".pdf", sep=""), "\n" )
  pdf(file=paste(out_dir,"MAplot_", comp, "_FDR_", fdr, ".pdf", sep=""), height=6, width=6)
  detags <- rownames(tTags$table) # n= has to be the number of significantly differentially expressed genes
  plotSmear(lrt,de.tags=detags, cex=0.5) # plot of fold change given CPM, red for those < FDR
  abline(h=c(-1,1), col="dodgerblue") # blue line at twofold change
  dev.off()
  # Merge annotations with DGE stats
  tTags_frame <- tTags$table # make data frame from tTags
  tTags_frame$transcript_id <- rownames(tTags_frame)
  join <- merge(tTags_frame, annot) # gets the intersection of detags and transcripts, i.e., the annotations for the DE transcripts
  # Some summary stats
  cat("percent all transcripts with ORFs:", sum(annot$prot_id != ".") * 100 / nrow(annot), "\n")
  cat("percent toptags with ORFs:", sum(join$prot_id != ".") * 100 / nrow(join), "\n")
  cat("percent all transcripts with Pfam or BlastX or BlastP annotation:", sum((annot$Pfam != ".") | (annot$Top_BLASTX_hit != ".") | (annot$Top_BLASTP_hit != ".")) * 100 / nrow(annot), "\n")
  cat("percent toptags with Pfam or BlastX or BlastP annotation:", sum((join$Pfam != ".") | (join$Top_BLASTP_hit != ".") | (join$Top_BLASTX_hit != ".")) * 100 / nrow(join), "\n")
  cat("number of genes upregulated, (FC > ", critFC, "):", sum(join$logFC > log2(critFC)), "\n")
  cat("number of genes downregulated, (FC < ", 1/critFC, "):", sum(join$logFC < -log2(critFC)), "\n")
  join <- join[ abs(join$logFC) > log2(critFC) , ]
  # write a csv and return the relevant dataframe
  if(onlyAnnot) {
    cat("saving .csv of up or downregulated genes, given critical FC, with only annotated sequences:", paste(out_dir, "topTags_", comp, "_", "critFC", critFC, "_", fdr, "FDR_only_w_Annot.csv", sep="" ))
    write.csv( join[ (join$Pfam != ".") | (join$Top_BLASTX_hit != ".") | (join$Top_BLASTP_hit != ".") , ], file=paste(out_dir, "topTags_", comp, "_", "critFC", critFC, "_", fdr, "FDR_only_w_Annot.csv", sep="" ) )
    return(join[ (join$Pfam != ".") | (join$Top_BLASTX_hit != ".") | (join$Top_BLASTP_hit != ".") , ])
  } else {
    cat("saving .csv of up/downregulated genes, given a critical FC:", paste(out_dir, "topTags_", comp, "_", "critFC", critFC, "_", fdr, "FDR.csv", sep="" ))
    write.csv( join, file=paste(out_dir, "topTags_", comp, "_", "critFC", critFC, "_", fdr, "FDR.csv", sep="" ) )
    return(join)    
  }
}


# display the coefficients
colnames(fit)

# get lists of differentially expressed genes from each comparison using getDEGs function
DEGlist_lrt_pH_1 <- get_DEGs(lrt_pH_1, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_pH_2 <- get_DEGs(lrt_pH_2, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_pH_3 <- get_DEGs(lrt_pH_3, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_temp_1 <- get_DEGs(lrt_temp_1, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_temp_2 <- get_DEGs(lrt_temp_2, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_temp_3 <- get_DEGs(lrt_temp_3, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_inter_0 <- get_DEGs(lrt_inter_0, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_inter_1 <- get_DEGs(lrt_inter_1, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_inter_2 <- get_DEGs(lrt_inter_2, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_inter_3 <- get_DEGs(lrt_inter_3, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_inter_4 <- get_DEGs(lrt_inter_4, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_inter_5 <- get_DEGs(lrt_inter_5, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_inter_6 <- get_DEGs(lrt_inter_6, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_inter_7 <- get_DEGs(lrt_inter_7, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_inter_8 <- get_DEGs(lrt_inter_8, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_inter_9 <- get_DEGs(lrt_inter_9, annot=transcripts, fdr=0.01, critFC=2)

RKC2_DEGlist_pH <- rbind(DEGlist_lrt_pH_1, DEGlist_lrt_pH_2, DEGlist_lrt_pH_3)
unique_pH_DE_transcripts <- unique(RKC2_DEGlist_pH$transcript_id)
length(unique_pH_DE_transcripts) # number of uniquely differentially expressed transcripts across pH comparisons
write.csv(unique_pH_DE_transcripts, paste(out_dir, "pH_comparison_DEG_list.csv", sep=""))

RKC2_DEGlist_temp <- rbind(DEGlist_lrt_temp_1, DEGlist_lrt_temp_2, DEGlist_lrt_temp_3)
unique_temp_DE_transcripts <- unique(RKC2_DEGlist_temp$transcript_id)
length(unique_temp_DE_transcripts) # number of uniquely differentially expressed transcripts across temp comparisons
write.csv(unique_temp_DE_transcripts, paste(out_dir, "temp_comparison_DEG_list.csv", sep=""))

RKC2_DEGlist_interact <- rbind(DEGlist_lrt_inter_0,DEGlist_lrt_inter_1,DEGlist_lrt_inter_2,DEGlist_lrt_inter_3,DEGlist_lrt_inter_4,DEGlist_lrt_inter_5,DEGlist_lrt_inter_6,DEGlist_lrt_inter_7,DEGlist_lrt_inter_8,DEGlist_lrt_inter_9)
unique_interact_DE_transcripts <- unique(RKC2_DEGlist_interact$transcript_id)
length(unique_interact_DE_transcripts) # number of uniquely differentially expressed transcripts across temp comparisons
write.csv(unique_interact_DE_transcripts, paste(out_dir, "interact_comparison_DEG_list.csv", sep=""))

RKC2_DEGlist_all <- rbind(DEGlist_lrt_temp_1, DEGlist_lrt_temp_2, DEGlist_lrt_temp_3, DEGlist_lrt_pH_1, DEGlist_lrt_pH_2, DEGlist_lrt_pH_3, DEGlist_lrt_inter_0,DEGlist_lrt_inter_1,DEGlist_lrt_inter_2,DEGlist_lrt_inter_3,DEGlist_lrt_inter_4,DEGlist_lrt_inter_5,DEGlist_lrt_inter_6,DEGlist_lrt_inter_7,DEGlist_lrt_inter_8,DEGlist_lrt_inter_9)
unique_all_DE_transcripts <- unique(RKC2_DEGlist_all$transcript_id)
length(unique_all_DE_transcripts) # number of uniquely differentially expressed transcripts across temp comparisons
write.csv(unique_all_DE_transcripts, paste(out_dir, "all_comparison_DEG_list.csv", sep=""))

##
#
# Get all the transcripts that were being used (for annotation of less than all the data)
#
##

ALL_RKC_transcripts = rownames(lrt_inter_0$table)
write.table(ALL_RKC_transcripts,"All_RKC_transcriptID.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)


#####
#
# write cluster files
#
#####


write.cluster.file <- function(all_fpkm, unique_DE_transcripts, outfile, cnames) {
  # make matrix of fpkm values using only those transcripts found in all targetList
  FPKM_all <- data.matrix(all_fpkm[(rownames(all_fpkm) %in% unique_DE_transcripts),])
  # log2 transform, mean center rows
  FPKM_all = log2(FPKM_all+1)
  centered_data = t(scale(t(FPKM_all), scale=F)) # center rows, mean substracted
  # make data frame of centered_data matrix
  centeredDF <- as.data.frame(centered_data)
  # concatenate column names from groupings to make the column names more informative
  colnames(centeredDF) <- cnames
  # add CloneID vector to centeredDF
  centeredDF$CloneID <- rownames(centered_data)
  # get vector of cloneID
  cloneIDs <- centeredDF$CloneID
  ### generate a new vector of topBlastHits for those cloneIDs
  # generate a subset data frame of the matches to transcript_id in the annotation dataframe; this makes the search in the loop below much faster
  annot_matches <- transcripts[(transcripts$transcript_id %in% cloneIDs),]
  # initialize an empty vector
  blasthits <- c()
  # iterate over cloneID in centeredDF, appending each topBlastHit to the blasthits vector
  for (cloneID in centeredDF$CloneID) {
    match <- annot_matches[annot_matches$transcript_id %in% cloneID,]
    blasthits <- c(blasthits, match$Top_BLASTX_hit[1])
  }
  # add blasthits vector to centeredDF "NAME" column
  centeredDF$NAME <- blasthits
  # reorder columns
  centeredDF <- centeredDF[,c(ncol(centeredDF)-1,ncol(centeredDF),1:(ncol(centeredDF)-2))] # more generalized form
  # save centeredDF 
  write.table(centeredDF, file=outfile, quote=FALSE, row.names=FALSE, sep="\t")
}



unique_pH_DE_transcripts
unique_temp_DE_transcripts
unique_all_DE_transcripts
unique_interact_DE_transcripts


# concatenate column names from sample ID and groupings to make the column names more informative
# function to concatenate names
concat_name <- function(x){
  paste(x[1],x[2],sep="_")
}

# make a data frame of the sample ID and groupings
sampleDF <- data.frame(sample=colnames(all_fpkm), treatment=group)
sampleNames <- apply(sampleDF, 1, concat_name)
sampleNames

write.cluster.file(all_fpkm, unique_pH_DE_transcripts, paste(out_dir, "centered_data_pH_comparison.txt", sep=""), sampleNames)
write.cluster.file(all_fpkm, unique_temp_DE_transcripts, paste(out_dir, "centered_data_temp_comparison.txt", sep=""), sampleNames)
write.cluster.file(all_fpkm, unique_all_DE_transcripts, paste(out_dir, "centered_data_all_comparison.txt", sep=""), sampleNames)
write.cluster.file(all_fpkm, unique_interact_DE_transcripts, paste(out_dir, "centered_data_interact_comparison.txt", sep=""), sampleNames)


# # sort by logFC
# sorted_tags <- annotated_DEGs_pHlow[with(annotated_DEGs_pHlow, order(logFC)), ]
# 
# # write a .csv based on a threshold FC, in this case, ( FC < 0.5 ) and ( FC > 2 )
# write.csv(annotated_DEGs_pHlow[ abs(annotated_DEGs_pHlow$logFC) > 1 , ], file=paste( out_dir, "topTags_DEGlist_pHlow_FDR0_00001_PfamHits_logFC_greater_than_2.csv", sep="" )
# 


#######
#
# Read in DE gene table to redo MDS on just the differentially expressed genes
#
#####

newMDSdata=read.delim(file.choose(),sep="\t")
plotMDS(newMDSdata[,c(3:47)] , main = "MDS Plot of DE genes", labels = group, cex=0.7)




##########
#
# Functions for genewise data analysis
#
##########



#####
#
# JS plyr version
#
#####
pH2=pH # use as.factor from group file, or make this as a factor
Temperature # use as.factor from group file, or make this as a factor

library(plyr)
library(ggplot2)

gene="comp102584_c11_seq1"

get_gene_data <- function(gene, all_fpkm) {
  gene_data <- data.frame(fpkm=as.numeric(all_fpkm[gene,]), pH=pH2, Temperature=temp)
  fpkm_mean<-ddply(gene_data,.(pH2,Temperature), function(d) mean(d$fpkm)) # if group file has columns called pH and Temperature
  SD<-ddply(gene_data,.(pH2,Temperature), function(d) sd(d$fpkm)) # if group file has columns called pH and Temperature
  names(fpkm_mean)=c("pH","Temperature","Mean") # rename
  names(SD)=c("pH","Temperature","SD") # rename
  out_frame=merge(fpkm_mean,SD) #make one data frame
  out_frame=out_frame[,c(3,4,1,2)] # reorder columns
  return(out_frame)
}

gene_data <- get_gene_data("comp102584_c11_seq1", all_fpkm)  #the gene in all_fpkm for which data will be returned

grep("comp102584_c11_seq1", all_fpkm)



####
# Function: get_gene_data
# examine expression of an individual gene based on transcript_id
# Scott's method
####
# 
# 
# 
# get_gene_data <- function(gene, fpkm) {
#   A_mean <- apply(fpkm[gene,1:5], 1, mean)
#   B_mean <- apply(fpkm[gene,6:10], 1, mean)
#   C_mean <- apply(fpkm[gene,11:15], 1, mean)
#   D_mean <- apply(fpkm[gene,16:20], 1, mean)
#   E_mean <- apply(fpkm[gene,21:25], 1, mean)
#   F_mean <- apply(fpkm[gene,26:30], 1, mean)
#   G_mean <- apply(fpkm[gene,31:35], 1, mean)
#   H_mean <- apply(fpkm[gene,36:40], 1, mean)
#   J_mean <- apply(fpkm[gene,41:45], 1, mean)
#   A_sd <- apply(fpkm[gene,1:5], 1, sd)
#   B_sd <- apply(fpkm[gene,6:10], 1, sd)
#   C_sd <- apply(fpkm[gene,11:15], 1, sd)
#   D_sd <- apply(fpkm[gene,16:20], 1, sd)
#   E_sd <- apply(fpkm[gene,21:25], 1, sd)
#   F_sd <- apply(fpkm[gene,26:30], 1, sd)
#   G_sd <- apply(fpkm[gene,31:35], 1, sd)
#   H_sd <- apply(fpkm[gene,36:40], 1, sd)
#   J_sd <- apply(fpkm[gene,41:45], 1, sd)
#   fpkm_mean <- c(A_mean, B_mean, C_mean, D_mean, E_mean, F_mean, G_mean, H_mean, J_mean)
#   SD <- c(A_sd, B_sd, C_sd, D_sd, E_sd, F_sd, G_sd, H_sd, J_sd)
#   names <- c("A", "B", "C", "D", "E", "F", "G", "H", "J")
#   pH_grp <- c("amb", "low", "amb", "low", "mid", "amb", "mid", "low", "mid")
#   temp_grp <- c("13C", "11C", "11C", "14C", "13C", "14C", "14C", "13C", "11C")
#   out_frame <- data.frame(fpkm_mean, SD, pH_grp, temp_grp, names)
#   out_frame$pH_grp <- factor( out_frame$pH_grp, levels = c( "amb", "mid", "low" ) )
#   return(out_frame)
# }

# get fpkm data for a specific gene
#gene_data <- get_gene_data("comp106087_c0_seq1", all_fpkm)


#####
# Function: get_annot
# gets annotation for a particular gene
#####

get_annot <- function(gene_in, annotation) {
  L <- annotation$transcript_id == gene_in
  cat("transcript: ", as.vector(annotation[L,]$transcript_id), "\n")
  cat("top blast hit: ", as.vector(annotation[L,]$Top_BLASTX_hit ), "\n")
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

# plot_genewise_fpkm <- function(gene_in, plot_title, add_title) {
#   k <- ggplot(gene_in, aes( x = temp_grp, y = fpkm_mean, ymax = fpkm_mean + SD, ymin = fpkm_mean - SD, fill = temp_grp ))
#   k + geom_bar(stat="identity") + facet_grid( . ~ pH_grp ) + geom_errorbar( width = 0.3) + scale_fill_brewer(type="qual", palette="Paired") + labs(title=paste(plot_title, add_title))
#   k + geom_bar(stat="identity") + facet_grid( . ~ pH_grp ) + geom_errorbar( width = 0.3) + scale_fill_brewer(type="qual", palette="Paired") + labs(title=paste(plot_title, add_title))
#   ggsave(filename = paste(out_dir, plot_title, "_", add_title, "_genewise_FPKM", ".pdf", sep =""), width = 6, height = 6)
#   }

# test plotting fpkm data
#plot_genewise_fpkm(gene_data, "titleA", "titleB")

#####
# function: plot_interaction
# makes an interaction plot for a given gene
#####

plot_interaction <- function(gene_data, plot_title, add_title) {
  h <- ggplot(gene_data, aes(x = pH, y = Mean, ymax = Mean + SD, ymin = Mean - SD, color=Temperature, fill=Temperature))
  h + geom_errorbar(color="black", width=0.1) + geom_point(aes(shape=Temperature, fill=Temperature),size = 4) + geom_line(aes(group=Temperature), size=0.5, linetype=1) + labs(title=paste(plot_title, add_title), x = "pH treatment", y = "mean FPKM") + scale_color_manual(values=c("#3399FF", "#66CC00", "#CC0066"))
  ggsave(filename = paste(out_dir, plot_title, "_", add_title, "_interactionPlot_FPKM", ".pdf", sep =""), width = 6, height = 6)
  }

#plot_interaction(gene_data,"TITLEA","TITLEB")
#####
# function: GOI_data
# master function to get data and plots from Gene Of Interest
#   gene_name: the "transcript_id" for your gene of interest
#   annotation: the full trinotate annotation report as a data frame
#   fpkm_table: the table of FPKM values 
#####

GOI_data <- function(gene_name, annotation, fpkm_table, graph_title) {
  gene_data <- get_gene_data(gene_name, fpkm_table)
  get_annot(gene_name, annotation)
  #plot_genewise_fpkm(gene_data, gene_name, graph_title)
  plot_interaction(gene_data, gene_name, graph_title)
}

# generate genewise plots based on a gene of interest
GOI_data("comp100107_c0_seq1", transcripts, all_fpkm, "Carbonic Anhydrase")

###
#
#  Example: read in a list of genes of interest, and make plots
#
###

desired_graphs = read.csv(file.choose())  #a text file with two columns, the geneID (gene_name) and a descriptor (graph_title)

outdir = "/Users/jonathonstillman/Documents/RKC_Transcriptomics/RKC2_Juveniles/GOIplots_temperature/"

for (i in 1:nrow(desired_graphs)) {
  GOI_data(desired_graphs[i,1], transcripts, all_fpkm, desired_graphs[i,2])
}



###########
#
# Functions for plotting clusters as interactions plots.
#
############

# read in the clustered median-centered log-2 FPKM data (from the .cdt file made by Cluster 3.0 software).  The file has had the E-weights removed, a column of clusters added, and saved as a .csv

clustereddata=read.csv(file.choose())

# prepare data for plotting using reshape
library(reshape)
clusterdata2=melt(clustereddata[,c(3:48)],id.var=c("Cluster"))
library(plyr)
library(ggplot2)
head(clusterdata2)
clusterdata2$group = substr(clusterdata2$variable,1,6)


Mean_by_cluster<-ddply(clusterdata2,.(Cluster,variable), function(d) mean(d$value)) # if group file has columns called pH and Temperature
sd_by_cluster<-ddply(clusterdata2,.(Cluster,variable), function(d) sd(d$value)) # if group file has columns called pH and Temperature
names(Mean_by_cluster)=c("cluster","group","Mean") # rename
names(sd_by_cluster)=c("cluster","group","SD") # rename
plotting_data_by_cluster=merge(Mean_by_cluster,sd_by_cluster) #make one data frame

#REDO WITH PROPER GLOBAL STANDARD DEVIATION AND STANDARD ERROR + 1.98*SE (95% conf. int.)

head(plotting_data_by_cluster)
pHbyCluster=c(rep(c(rep("8.0",5),rep("7.5",5),rep("8.0",5),rep("7.5",5),rep("7.8",5),rep("8.0",5),rep("7.8",5),rep("7.5",5),rep("7.8",5)),6))
TempbyCluster=c(rep(c(rep("13",5),rep("11",5),rep("11",5),rep("14",5),rep("13",5),rep("14",5),rep("14",5),rep("13",5),rep("11",5)),6))
Library_cluster_plot=c(rep(c(rep("A",5),rep("B",5),rep("C",5),rep("D",5),rep("E",5),rep("F",5),rep("G",5),rep("H",5),rep("J",5)),6))
new_plotting_data_by_cluster=cbind(plotting_data_by_cluster,pHbyCluster,TempbyCluster,Library_cluster_plot)
names(new_plotting_data_by_cluster)=c("Cluster","Library","Mean","SD","pH","Temperature","Group")
head(new_plotting_data_by_cluster)

# new plot data with data averaged per group at each pH and Temperature and cluster
new_plotting_data_by_cluster_avgd<-ddply(clusterdata2,.(Cluster,group), function(d) mean(d$value)) # if group file has columns called pH and Temperature
sd_by_cluster_avgd<-ddply(clusterdata2,.(Cluster,group), function(d) sd(d$value)) # if group file has columns called pH and Temperature
library(plotrix)
se_by_cluster_avgd<-ddply(clusterdata2,.(Cluster,group), function(d) std.error(d$value)) 
head(new_plotting_data_by_cluster_avgd)
names(new_plotting_data_by_cluster_avgd)=c("Cluster","Group","Mean") # rename
names(sd_by_cluster_avgd)=c("Cluster","Group","SD") # rename
names(se_by_cluster_avgd)=c("Cluster","Group","SE")
plotting_data_by_cluster_avgd=merge(new_plotting_data_by_cluster_avgd,sd_by_cluster_avgd)
plotting_data_by_cluster_avgd=merge(plotting_data_by_cluster_avgd,se_by_cluster_avgd)
plotting_data_by_cluster_avgd$pct=plotting_data_by_cluster_avgd$SE*1.98
#make one 
head(plotting_data_by_cluster_avgd)
plotting_data_by_cluster_avgd$pH = c(rep(c(rep("8.1",3),rep("7.5",3),rep("7.8",3))))
plotting_data_by_cluster_avgd$Temperature = c(rep(c("10.5","12.5","14.5"),3))
head(plotting_data_by_cluster_avgd)

#plot with SD
c <- ggplot(plotting_data_by_cluster_avgd, aes(x = pH, y = Mean, ymax = Mean + SD, ymin = Mean - SD, color=Temperature, fill=Temperature,group=Group))
c + geom_errorbar(color="black", width=0.1) + geom_point(aes(shape=Temperature, fill=Temperature),size = 4) + geom_line(aes(group=Temperature), size=0.5, linetype=1) + labs(title="Clusters (error bars are 1SD)",x = "pH treatment", y = "Median centered log2 FPKM") + scale_color_manual(values=c("#3399FF", "#66CC00", "#CC0066"))+facet_grid(~Cluster)

#plot with SEM
d <- ggplot(plotting_data_by_cluster_avgd, aes(x = pH, y = Mean, ymax = Mean + SE, ymin = Mean - SE, color=Temperature, fill=Temperature,group=Group))
d + geom_errorbar(color="black", width=0.1) + geom_point(aes(shape=Temperature, fill=Temperature),size = 4) + geom_line(aes(group=Temperature), size=0.5, linetype=1) + labs(title="Clusters (error bars are 1SEM)",x = "pH treatment", y = "Median centered log2 FPKM") + scale_color_manual(values=c("#3399FF", "#66CC00", "#CC0066"))+facet_grid(~Cluster)

#plot with 95% conf intervals
e <- ggplot(plotting_data_by_cluster_avgd, aes(x = pH, y = Mean, ymax = Mean + pct, ymin = Mean - pct, color=Temperature, fill=Temperature,group=Group))
e + geom_errorbar(color="black", width=0.1) + geom_point(aes(shape=Temperature, fill=Temperature),size = 3) + geom_line(aes(group=Temperature), size=0.5, linetype=1) + labs(title="Juvenile Cluster 7 (error bars are 95% conf. intervals)",x = "pH treatment", y = "Median centered log2 FPKM") + scale_color_manual(values=c("#3399FF", "#66CC00", "#CC0066"))+theme_bw()#+facet_grid(~Cluster)

## NOW - do a 2-way ANOVA (pH and T) for each Cluster with TUKEY HSD for PW comparisions
# subset the data by cluster
Cluster1_data=subset(clusterdata2,clusterdata2$Cluster==1)
Cluster1_data$pH = substr(Cluster1_data$group,1,3)
Cluster1_data$Temperature = substr(Cluster1_data$group,4,5)

Cluster2_data=subset(clusterdata2,clusterdata2$Cluster==2)
Cluster2_data$pH = substr(Cluster2_data$group,1,3)
Cluster2_data$Temperature = substr(Cluster2_data$group,4,5)

Cluster3_data=subset(clusterdata2,clusterdata2$Cluster==3)
Cluster3_data$pH = substr(Cluster3_data$group,1,3)
Cluster3_data$Temperature = substr(Cluster3_data$group,4,5)

Cluster4_data=subset(clusterdata2,clusterdata2$Cluster==4)
Cluster4_data$pH = substr(Cluster4_data$group,1,3)
Cluster4_data$Temperature = substr(Cluster4_data$group,4,5)

Cluster5_data=subset(clusterdata2,clusterdata2$Cluster==5)
Cluster5_data$pH = substr(Cluster5_data$group,1,3)
Cluster5_data$Temperature = substr(Cluster5_data$group,4,5)

Cluster6_data=subset(clusterdata2,clusterdata2$Cluster==6)
Cluster6_data$pH = substr(Cluster6_data$group,1,3)
Cluster6_data$Temperature = substr(Cluster6_data$group,4,5)

# perform ANOVA for each cluster to ID genes that differed by group.
Cluster1AOV=aov(value~pH*Temperature, data=Cluster1_data)
TukeyHSD(Cluster1AOV)
write.csv(HSD1$`pH:Temperature`,"TukeyHSD1.csv")
Cluster2AOV=aov(value~pH*Temperature, data=Cluster2_data)
HSD2= TukeyHSD(Cluster2AOV)
write.csv(HSD2$`pH:Temperature`,"TukeyHSD2.csv")
Cluster3AOV=aov(value~pH*Temperature, data=Cluster3_data)
HSD3= TukeyHSD(Cluster3AOV)
write.csv(HSD3$`pH:Temperature`,"TukeyHSD3.csv")
Cluster4AOV=aov(value~pH*Temperature, data=Cluster4_data)
HSD4= TukeyHSD(Cluster4AOV)
write.csv(HSD4$`pH:Temperature`,"TukeyHSD4.csv")
Cluster5AOV=aov(value~pH*Temperature, data=Cluster5_data)
HSD5= TukeyHSD(Cluster5AOV)
write.csv(HSD5$`pH:Temperature`,"TukeyHSD5.csv")
Cluster6AOV=aov(value~pH*Temperature, data=Cluster6_data)
HSD6= TukeyHSD(Cluster6AOV)
write.csv(HSD6$`pH:Temperature`,"TukeyHSD6.csv")


#######
#
# Make a table with specific gene text (e.g., "resilin")
#
#######

#subset the data to make a new table just containing text string.
resilindata=subset(newMDSdata, regexpr("resilin", newMDSdata$NAME) > 0)

cuticledata=subset(newMDSdata, regexpr("cuticle", newMDSdata$NAME) > 0)
library(reshape)
resilindatamelt=melt(resilindata,id.vars="CloneID",measure.vars=c(resilindata[3,47]))
resilindatamelt=melt(resilindata)
colnames(resilindata)
pHlong=c(rep("8.0",255),rep("7.5",255),rep("8.0",255),rep("7.5",255),rep("7.8",255),rep("8.0",255),rep("7.8",255),rep("7.5",255),rep("7.8",255))
Templong=c(rep("13",255),rep("11",255),rep("11",255),rep("14",255),rep("13",255),rep("14",255),rep("14",255),rep("13",255),rep("11",255))

resilindataplot=cbind(resilindatamelt,pHlong,Templong)
head(resilindataplot)
library(plyr)
library(ggplot2)
Mean<-ddply(resilindataplot,.(pHlong,Templong), function(d) mean(d$value)) # if group file has columns called pH and Temperature
sd<-ddply(resilindataplot,.(pHlong,Templong), function(d) sd(d$value)) # if group file has columns called pH and Temperature
names(Mean)=c("pH","Temperature","Mean") # rename
names(sd)=c("pH","Temperature","SD") # rename
plotting_data=merge(Mean,sd) #make one data frame
plotting_data$pH = c(rep(7.5,3),rep(7.8,3),rep(8.0,3))
a <- ggplot(plotting_data, aes(x = pH, y = Mean, ymax = Mean + SD, ymin = Mean - SD, color=Temperature, fill=Temperature))
a + geom_errorbar(color="black", width=0.02) + geom_point(aes(shape=Temperature, fill=Temperature),size = 4) + geom_line(aes(group=Temperature), size=0.5, linetype=1) + labs(title="all pro-resilins",x = "pH treatment", y = "mean FPKM") + scale_color_manual(values=c("#3399FF", "#66CC00", "#CC0066"))+scale_x_continuous(breaks=c(7.5,7.8,8.0))


## make an interaction plot of the unclustered cluster:
unclustdata=read.csv(file.choose())
# prepare data for plotting using reshape
library(reshape)
clusterdata2=melt(unclustdata[,c(3:48)],id.var=c("Cluster"))
library(plyr)
library(ggplot2)
clusterdata2$group = substr(clusterdata2$variable,1,6)
head(clusterdata2)

