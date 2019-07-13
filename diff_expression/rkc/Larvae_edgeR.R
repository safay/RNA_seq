#!/usr/bin/R

# Scott Fay 
# 2013 Oct 15
# updated 10 April 2014

# to install packages:
# install.packages("ggplot2")
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("edgeR")

library(ggplot2)
library(edgeR)

#####
#
# User-defined variables
#
#####

# set working directory; must only contain directories containing eXpress output, no other files or directories 
setwd("/Users/Shared/RKC_project/RKC3_RKC4/express_output/")
# set output directory
out_dir <- "/Users/Shared/RKC_project/RKC3_RKC4/DGE_analysis_output/"
# get annotation file
# "trinotate_annotation_report.txt" is a Trinotate output excel file, exported as tab-delimited
transcripts <- read.delim("/Users/scottfay/Dropbox/Red_King_Crab_Project/RKC_transcriptome_and_annotation/trinotate_annotation_report.txt",stringsAsFactors=FALSE)

# define groups
groups <- read.delim("../RKC3_RKC4_groups.txt", colClasses = c("character", rep("factor",4)))
str(groups)

##########
#
# Import eXpress expression estimates, total counts and FPKM values
#
##########

# creates a vector of the sequencing library directory names, "libs"
libs <- list.files()

# loads in the results.xprs files
# each is in a unique directory in the pwd
for (dirname in libs) {
  cat('reading in', dirname, '\n')
  tempframe <- read.delim(paste('./', dirname, '/results.xprs', sep=''))
  assign(dirname, tempframe[order(tempframe$target_id),]) # generate dataframe, organizing rows by target_id
}

# Make an FPKM dataframe for heatmap and plots, "all_fpkm"
# initialize column of transcript names and the first data column
all_fpkm <- data.frame(target_id=get(libs[1])$target_id, fpkm=get(libs[1])$fpkm)
# change the column name
colnames(all_fpkm)[colnames(all_fpkm)=="fpkm"] <- paste(libs[1],"_fpkm",sep="") 
# now iterate over the remaining columns...
for (name in libs[2:length(libs)]) {
  # makes a new column
  all_fpkm$fpkm <- get(name)$fpkm
  # renames that new column
  colnames(all_fpkm)[colnames(all_fpkm)=="fpkm"] <- paste(name,"_fpkm",sep="")
}

# Make a counts dataframe for EdgeR analysis, "all_count"
# initialize column of transcript names and the first data column
all_count <- data.frame(target_id=get(libs[1])$target_id, tot_counts=get(libs[1])$tot_counts)
# change the column name
colnames(all_count)[colnames(all_count)=="tot_counts"] <- paste(libs[1],"_count",sep="")
# now iterate over the remaining columns...
for (name in libs[2:length(libs)]) {
  # makes a new column
  all_count$tot_counts <- get(name)$tot_counts
  # renames that new column
  colnames(all_count)[colnames(all_count)=="tot_counts"] <- paste(name,"_count",sep="")
}

# clean up environment
rm(tempframe, name, dirname)
rm(list = libs) 

# create row names from target_id
rownames(all_count) <- all_count$target_id
rownames(all_fpkm) <- all_fpkm$target_id
# remove target_id column
all_count$target_id <- NULL
all_fpkm$target_id <- NULL

# Rename columns
colnames(all_count)
colnames(all_count) <- c("1", "10", "11", "112", "123", "13", "134", "145", "156", "167", "178", "189", "200", "211", "222", "233", "24", "244", "255", "266", "277", "288", "299", "3", "4", "46", "5", "57", "6", "79", "9", "90")
colnames(all_fpkm) <- c("1", "10", "11", "112", "123", "13", "134", "145", "156", "167", "178", "189", "200", "211", "222", "233", "24", "244", "255", "266", "277", "288", "299", "3", "4", "46", "5", "57", "6", "79", "9", "90")
# now the data is ready

# remove spike option
# retreive the spike IDs generated using blat
# something like this:
# blat /Users/scottfay/Google\ Drive/RNA_seq/Spikes/spikes_1_8.fa ~/annotation/Dicosmo/Dicosmo.Trinity.fasta ~/annotation/Dicosmo/spikes_blat_output
# and then just manually edit out the spike transcript IDs into a separate file
#spikes <- scan("/Users/scottfay/annotation/Dicosmo/Dicosmo_spike_transcript_IDs.txt", what = "character")
# remove the spikes transcripts from the tables
#all_count <- all_count[!rownames(all_count) %in% spikes, ]
#all_fpkm <- all_fpkm[!rownames(all_fpkm) %in% spikes, ]

### as a test, remove samples B1 and B4...
#new_count <- all_count[,c(1:5,7,8,10:20)]
#all_count <- new_count
#head(all_count)


#####
#
# filter out transcripts by mean count per gene
# NOTE: another kind of filtering, that takes into account library size, is used below.  Prefer that method over this one.
#
####

filterVal <- 0    # filters out mean count per gene of less than this value
means <- rowMeans(all_count)
filter <- means >= filterVal
table(filter)
filtered.all.count <- all_count[filter,]
dim(filtered.all.count)

#####
#
# Plot total number of mapped reads to known genes: look for outliers of library size
#
#####

library(RColorBrewer)
colors <- brewer.pal(9, "Set1")
boxplot(log2(all_count+1), las=2) # before filter
boxplot(log2(filtered.all.count+1), las=2) # after filter

#####
#
# Drop any poor quality libraries
#
#####

drops <- c("211","233")  # list the library names to be removed here
filtered.all.count <- filtered.all.count[,!(names(filtered.all.count) %in% drops)]
all_fpkm <- all_fpkm[,!(names(all_fpkm) %in% drops)]
groups <- groups[!(groups$Sample.ID %in% drops),]
boxplot(log2(filtered.all.count+1), las=2) # visualize again

# make three sets of libraries for downstream analysis:

# 1. newly hatched larvae
new.hatch.groups <- groups[1:8,]
new.hatch.libs <- filtered.all.count[ , names(filtered.all.count) %in% new.hatch.groups$Sample.ID]
new.hatch.libs <- new.hatch.libs[,new.hatch.groups$Sample.ID]  # reorder so groups and libs have the same order

# 2. just 7-day larvae
seven.day.groups <- groups[9:30,]
seven.day.libs <- filtered.all.count[ , names(filtered.all.count) %in% seven.day.groups$Sample.ID]
seven.day.libs <- seven.day.libs[,seven.day.groups$Sample.ID]  # reorder libs columns so groups and libs have the same order

# 3. both together
all.libs <- filtered.all.count
all.groups <- groups
all.libs <- all.libs[,all.groups$Sample.ID]  # reorder libs columns so groups and libs have the same order

#####
#
# Make EdgeR objects and normalize
#
#####

x <- DGEList(counts=new.hatch.libs)
y <- DGEList(counts=seven.day.libs) 
z <- DGEList(counts=all.libs)  

# TMM normalization
x.tmm <- calcNormFactors(x, method="TMM")
y.tmm <- calcNormFactors(y, method="TMM")
z.tmm <- calcNormFactors(z, method="TMM")

# upper quartile normalization
x.uq <- calcNormFactors(x, method="upperquartile")
y.uq <- calcNormFactors(y, method="upperquartile")
z.uq <- calcNormFactors(z, method="upperquartile")

#####
#
# MDS plots
#
#####

plotMDS(x , main = "MDS Plot for Count Data", labels = colnames( x$counts ), cex=0.7)
plotMDS(y , main = "MDS Plot for Count Data", labels = colnames( y$counts ), cex=0.7)
plotMDS(z , main = "MDS Plot for Count Data", labels = colnames( z$counts ), cex=0.7)

plotMDS(x.tmm , main = "MDS Plot for Count Data", labels = colnames( x$counts ), cex=0.7)
plotMDS(y.tmm , main = "MDS Plot for Count Data", labels = colnames( y$counts ), cex=0.7)
plotMDS(z.tmm , main = "MDS Plot for Count Data", labels = colnames( z$counts ), cex=0.7)

plotMDS(x.uq , main = "MDS Plot for Count Data", labels = colnames( x$counts ), cex=0.7)
plotMDS(y.uq , main = "MDS Plot for Count Data", labels = colnames( y$counts ), cex=0.7)
plotMDS(z.uq , main = "MDS Plot for Count Data", labels = colnames( z$counts ), cex=0.7)

# other labels:
plotMDS(y.uq , main = "MDS Plot for Count Data", labels = colnames( y$counts ), cex=0.7)
plotMDS(y.uq , main = "MDS Plot for Count Data", labels = y$samples$group, cex=0.7)
plotMDS(y.uq , main = "MDS Plot for Count Data", labels = groups[9:30,]$mother, cex=0.7)

# more dimensions, defined in dim.plot()
plotMDS(y.uq , main = "MDS Plot for Count Data", labels = y$samples$group, , dim.plot=c(3,4), cex=0.7)

#####
# 
# Test of differential expression
# Analysis of the seven day larvae, in a blocked design to see effect of pH on the larvae while factoring out maternal effect
#
#####

# define the statistical model design matrix using the model.matrix function
mom <- factor(groups[9:30,]$mother)
mom <- factor(mom, levels=c("1594", "1572", "1603"))
pH <- factor(groups[9:30,]$larval.pH)
pH <- factor(pH, levels=c("8", "7.8", "7.5")) # get the factors in the "right" order
# this follows the blocking model in the EdgeR documentation 3.4.2, pg 32.
design <- model.matrix(~mom+pH, data=y.uq$samples)  #
design

# FILTERING: filter out transcripts that do not have at least five reads out of 1,000,000 reads mapped in at least 3 samples
# NOTE: BCV values look much better when you do this kind of filtering than just mean > a certain threshhold.
y.uq <- y.uq[rowSums(1e+06 * y.uq$counts/expandAsMatrix(y.uq$samples$lib.size, dim(y.uq)) > 5) >= 3, ]

# reset the library sizes after filtering
y.uq$samples$lib.size <- colSums(y.uq$counts)

# estimate dispersion for GLM fit, common, trended, and tagwise
# see edgeR docs for references on dispersion estimates, e.g., "?estimateGLMCommonDisp"
y.uq <- estimateGLMCommonDisp(y.uq,design) 
# estimate trended dispersion for use in tagwise dispersion estimate
y.uq <- estimateGLMTrendedDisp(y.uq,design)
# estimate tagwise dispersion to be used in glmFit()
y.uq <- estimateGLMTagwiseDisp(y.uq,design)
# plot genewise biological coefficient of variation against gene abundance
plotBCV(y.uq)
# make plot as a pdf
pdf(file=paste(out_dir, "Day7_BCV_plot.pdf", sep=""), height=6, width=6)
plotBCV(y.uq, main = "Day7_Biological Coefficient of Variation")
dev.off()
# MDS plot again
plotMDS(y.uq , main = "MDS Plot for Count Data", labels = colnames( y.uq$counts ), cex=0.7)
# make a pdf
#pdf(file=paste(out_dir, "MDS_plot.pdf", sep=""), height=6, width=6)
#plotMDS(y.uq , main = "MDS Plot for Count Data", labels = colnames( y$counts ), cex=0.7)
#dev.off()

# GLM
# fit the negative binomial generalized linear model (GLM) for each tag, creating a new fit object
fit <- glmFit(y.uq, design)
colnames(fit)

#####
# 
# LRT to show the top genes
#
#####

## pH comparisons
lrt_pH_1 <- glmLRT(fit, coef=4)
lrt_pH_2 <- glmLRT(fit, coef=5)
lrt_pH_3 <- glmLRT(fit, contrast=c(0,0,0,-1,1))

## mother comparisons
lrt_mom_1 <- glmLRT(fit, coef=2)
lrt_mom_2 <- glmLRT(fit, coef=3)
lrt_mom_3 <- glmLRT(fit, contrast=c(0,-1,1,0,0))


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
DEG <- summary(decideTestsDGE(lrt_pH_1, p=0.05, adjust="BH"))

get_DEGs <- function(lrt, annot, fdr=0.05, critFC=2, onlyAnnot=FALSE) {
  DEG <- summary(decideTestsDGE(lrt, p=fdr, adjust="BH")) # gives numbers of genes up and downregulated at FDR < pval
  Num_DEG <- (DEG[1] + DEG[-1])[2] # gets total number of DEGs based on the FDR
  if(Num_DEG==0){
    cat("There are no statistically significant differences in gene expression")
    break
  }
  tTags <- topTags(lrt, n=Num_DEG) # gets a list of DEGs, "topTags"
  cp <- unlist(strsplit(tTags$comparison, "[* ]"))
  comp <- paste(cp[2], "v", cp[4], sep="") # this and the preceding line make a cleaner text format for the "comparison," used in filenames.  This gets  screwed up when you use the glmLRT(fit, coef=4) style instead of the glmLRT(fit, contrast=c(0,0,0,-1,1)) style
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

# get lists of differentially expressed genes among mothers from each comparison using getDEGs function
DEGlist_lrt_mom_1 <- get_DEGs(lrt_mom_1, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_mom_2 <- get_DEGs(lrt_mom_2, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt_mom_3 <- get_DEGs(lrt_mom_3, annot=transcripts, fdr=0.01, critFC=2)

RKC4_DEGlist_mom <- rbind(DEGlist_lrt_mom_1, DEGlist_lrt_mom_2, DEGlist_lrt_mom_3)
unique_mom_DE_transcripts <- unique(RKC4_DEGlist_mom$transcript_id)
length(unique_mom_DE_transcripts) # number of uniquely differentially expressed transcripts across mom comparisons
write.csv(unique_mom_DE_transcripts, paste(out_dir, "mom_comparison_DEG_list.csv", sep=""))

# get lists of differentially expressed genes across pH levels from each comparison using getDEGs function
DEGlist_lrt_pH_1 <- get_DEGs(lrt_pH_1, annot=transcripts, fdr=0.1, critFC=0) # 0 DE genes
DEGlist_lrt_pH_2 <- get_DEGs(lrt_pH_2, annot=transcripts, fdr=0.1, critFC=1)  # 2 genes!
DEGlist_lrt_pH_3 <- get_DEGs(lrt_pH_3, annot=transcripts, fdr=0.1, critFC=1) # 0 DE genes
write.csv(DEGlist_lrt_pH_2, paste(out_dir, "pH_comparison_DEG_list.csv", sep=""))


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
    blasthits <- c(blasthits, match$TopBlastHit[1])
  }
  # add blasthits vector to centeredDF "NAME" column
  centeredDF$NAME <- blasthits
  # reorder columns
  centeredDF <- centeredDF[,c(ncol(centeredDF)-1,ncol(centeredDF),1:(ncol(centeredDF)-2))] # more generalized form
  # save centeredDF 
  write.table(centeredDF, file=outfile, quote=FALSE, row.names=FALSE, sep="\t")
}





head(unique_mom_DE_transcripts)



# concatenate column names from sample ID and groupings to make the column names more informative
# function to concatenate names
concat_name <- function(x){
  paste(x[1],x[3],x[5],sep="_")
}

# make a data frame of the sample ID and groupings

sampleNames <- apply(groups[-c(1:8),], 1, concat_name)
sampleNames

write.cluster.file(all_fpkm[,c(seven.day.groups$Sample.ID)], unique_mom_DE_transcripts, paste(out_dir, "Day7_centered_data_mom_comparison.txt", sep=""), sampleNames)


# # sort by logFC
# sorted_tags <- annotated_DEGs_pHlow[with(annotated_DEGs_pHlow, order(logFC)), ]
# 
# # write a .csv based on a threshold FC, in this case, ( FC < 0.5 ) and ( FC > 2 )
# write.csv(annotated_DEGs_pHlow[ abs(annotated_DEGs_pHlow$logFC) > 1 , ], file=paste( out_dir, "topTags_DEGlist_pHlow_FDR0_00001_PfamHits_logFC_greater_than_2.csv", sep="" )
# 

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
pH # use as.factor from group file, or make this as a factor
mom # use as.factor from group file, or make this as a factor

library(plyr)
library(ggplot2)

gene="comp0_c0_seq1"

get_gene_data <- function(gene, all_fpkm) {
  gene_data <- data.frame(fpkm=as.numeric(all_fpkm[gene,]), pH=pH, Mom=mom)
  fpkm_mean<-ddply(gene_data,.(pH,mom), function(d) mean(d$fpkm)) # if group file has columns called pH and Mom
  SD<-ddply(gene_data,.(pH,mom), function(d) sd(d$fpkm)) # if group file has columns called pH and Mom
  names(fpkm_mean)=c("pH","Mom","Mean") # rename
  names(SD)=c("pH","Mom","SD") # rename
  out_frame=merge(fpkm_mean,SD) #make one data frame
  out_frame=out_frame[,c(3,4,1,2)] # reorder columns
  return(out_frame)
}

gene_data <- get_gene_data("comp0_c0_seq1", all_fpkm[,c(seven.day.groups$Sample.ID)])  #the gene in all_fpkm for which data will be returned

gene_data



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
  h <- ggplot(gene_data, aes(x = pH, y = Mean, ymax = Mean + SD, ymin = Mean - SD, color=Mom))
  h + geom_errorbar(color="black", width=0.1) + geom_point(aes(shape=Mom), color="black",size = 6) + geom_line(aes(group=Mom), size=0.5, linetype=1) + labs(title=paste(plot_title, add_title), x = "pH treatment", y = "mean FPKM") + scale_color_manual(values=c("#3399FF", "#66CC00", "#CC0066"))
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
GOI_data("comp100107_c0_seq1", transcripts, all_fpkm[,c(seven.day.groups$Sample.ID)], "Carbonic Anhydrase")

###
#
#  Example: read in a list of genes of interest, and make plots
#
###

desired_graphs = read.csv(file.choose())  #a text file with two columns, the geneID (gene_name) and a descriptor (graph_title)

for (i in 1:nrow(desired_graphs)) {
  GOI_data(desired_graphs[i,1], transcripts, all_fpkm[,c(seven.day.groups$Sample.ID)], desired_graphs[i,2])
}

#################################################################################################
#################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
#####
# 
# Test of differential expression
# Analysis of the Day 0 larvae, in a blocked design to see effect of pH on the larvae while factoring out maternal effect
#
#####

# define the statistical model design matrix using the model.matrix function
pH <- factor(new.hatch.groups$mother.pH)
pH <- factor(pH, levels=c("8", "7.8", "7.5")) # get the factors in the "right" order
# this follows the blocking model in the EdgeR documentation 3.4.2, pg 32.
design <- model.matrix(~pH, data=x.uq$samples)  #
design

# FILTERING: filter out transcripts that do not have at least five reads out of 1,000,000 reads mapped in at least 3 samples
# NOTE: BCV values look much better when you do this kind of filtering than just mean > a certain threshhold.
x.uq <- x.uq[rowSums(1e+06 * x.uq$counts/expandAsMatrix(x.uq$samples$lib.size, dim(x.uq)) > 5) >= 3, ]

# reset the library sizes after filtering
x.uq$samples$lib.size <- colSums(x.uq$counts)

# estimate dispersion for GLM fit, common, trended, and tagwise
# see edgeR docs for references on dispersion estimates, e.g., "?estimateGLMCommonDisp"
x.uq <- estimateGLMCommonDisp(x.uq,design) 
# estimate trended dispersion for use in tagwise dispersion estimate
x.uq <- estimateGLMTrendedDisp(x.uq,design)
# estimate tagwise dispersion to be used in glmFit()
x.uq <- estimateGLMTagwiseDisp(x.uq,design)
# plot genewise biological coefficient of variation against gene abundance
plotBCV(x.uq)
# make plot as a pdf
pdf(file=paste(out_dir, "Day0_BCV_plot.pdf", sep=""), height=6, width=6)
plotBCV(x.uq, main = "Biological Coefficient of Variation")
dev.off()
# MDS plot again
plotMDS(x.uq , main = "MDS Plot for Count Data", labels = colnames( x.uq$counts ), cex=0.7)
# make a pdf
#pdf(file=paste(out_dir, "MDS_plot.pdf", sep=""), height=6, width=6)
#plotMDS(x.uq , main = "MDS Plot for Count Data", labels = colnames( y$counts ), cex=0.7)
#dev.off()

# GLM
# fit the negative binomial generalized linear model (GLM) for each tag, creating a new fit object
fit <- glmFit(x.uq, design)
colnames(fit)

#####
# 
# LRT to show the top genes
#
#####

## pH comparisons
lrt_pH_1 <- glmLRT(fit, coef=2)
lrt_pH_2 <- glmLRT(fit, coef=3)
lrt_pH_3 <- glmLRT(fit, contrast=c(0,-1,1))


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
DEG <- summary(decideTestsDGE(lrt_pH_1, p=0.05, adjust="BH"))

get_DEGs <- function(lrt, annot, fdr=0.05, critFC=2, onlyAnnot=FALSE) {
  DEG <- summary(decideTestsDGE(lrt, p=fdr, adjust="BH")) # gives numbers of genes up and downregulated at FDR < pval
  Num_DEG <- (DEG[1] + DEG[-1])[2] # gets total number of DEGs based on the FDR
  if(Num_DEG==0){
    cat("There are no statistically significant differences in gene expression")
    break
  }
  tTags <- topTags(lrt, n=Num_DEG) # gets a list of DEGs, "topTags"
  cp <- unlist(strsplit(tTags$comparison, "[* ]"))
  comp <- paste(cp[2], "v", cp[4], sep="") # this and the preceding line make a cleaner text format for the "comparison," used in filenames.  This gets  screwed up when you use the glmLRT(fit, coef=4) style instead of the glmLRT(fit, contrast=c(0,0,0,-1,1)) style
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

# get lists of differentially expressed genes among mothers from each comparison using getDEGs function
DEGlist_lrt_pH_1 <- get_DEGs(lrt_pH_1, annot=transcripts, fdr=0.01, critFC=2) # no differences
DEGlist_lrt_pH_2 <- get_DEGs(lrt_pH_2, annot=transcripts, fdr=0.01, critFC=2) # no differences
DEGlist_lrt_pH_3 <- get_DEGs(lrt_pH_3, annot=transcripts, fdr=0.01, critFC=2) # 46 differences

#RKC4_DEGlist_pH <- rbind(DEGlist_lrt_pH_1, DEGlist_lrt_pH_2, DEGlist_lrt_pH_3)
#unique_pH_DE_transcripts <- unique(RKC4_DEGlist_pH$transcript_id)
#length(unique_pH_DE_transcripts) # number of uniquely differentially expressed transcripts across pH comparisons
write.csv(DEGlist_lrt_pH_3, paste(out_dir, "Day0_pH_comparison_DEG_list.csv", sep=""))

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
    blasthits <- c(blasthits, match$TopBlastHit[1])
  }
  # add blasthits vector to centeredDF "NAME" column
  centeredDF$NAME <- blasthits
  # reorder columns
  centeredDF <- centeredDF[,c(ncol(centeredDF)-1,ncol(centeredDF),1:(ncol(centeredDF)-2))] # more generalized form
  # save centeredDF 
  write.table(centeredDF, file=outfile, quote=FALSE, row.names=FALSE, sep="\t")
}

head(DEGlist_lrt_pH_3)

# concatenate column names from sample ID and groupings to make the column names more informative
# function to concatenate names
concat_name <- function(x){
  paste(x[1],x[2],x[5],sep="_")
}

# make a data frame of the sample ID and groupings

sampleNames <- apply(new.hatch.groups, 1, concat_name)
sampleNames

larvae_0day_DE_transcripts <- DEGlist_lrt_pH_3$transcript_id

write.cluster.file(all_fpkm[,c(new.hatch.groups$Sample.ID)], larvae_0day_DE_transcripts, paste(out_dir, "Day0_centered_data_pH_comparison.txt", sep=""), sampleNames)


# # sort by logFC
# sorted_tags <- annotated_DEGs_pHlow[with(annotated_DEGs_pHlow, order(logFC)), ]
# 
# # write a .csv based on a threshold FC, in this case, ( FC < 0.5 ) and ( FC > 2 )
# write.csv(annotated_DEGs_pHlow[ abs(annotated_DEGs_pHlow$logFC) > 1 , ], file=paste( out_dir, "topTags_DEGlist_pHlow_FDR0_00001_PfamHits_logFC_greater_than_2.csv", sep="" )
# 

##########
#
# Functions for genewise data analysis: Works!
#
##########



#####
#
# JS plyr version
#
#####
pH # use as.factor from group file, or make this as a factor
mom # use as.factor from group file, or make this as a factor

library(plyr)
library(ggplot2)

gene="comp0_c0_seq1"

get_gene_data <- function(gene, all_fpkm) {
  gene_data <- data.frame(fpkm=as.numeric(all_fpkm[gene,]), pH=pH)
  fpkm_mean<-ddply(gene_data,.(pH), function(d) mean(d$fpkm)) # if group file has columns called pH 
  SD<-ddply(gene_data,.(pH), function(d) sd(d$fpkm)) # if group file has columns called pH 
  names(fpkm_mean)=c("pH","Mean") # rename
  names(SD)=c("pH","SD") # rename
  out_frame=merge(fpkm_mean,SD) #make one data frame
  out_frame=out_frame[,c(2,3,1)] # reorder columns
  return(out_frame)
}

gene_data <- get_gene_data("comp0_c0_seq1", all_fpkm[,c(new.hatch.groups$Sample.ID)])  #the gene in all_fpkm for which data will be returned

gene_data



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
  k <- ggplot(gene_in, aes( x = pH, y = Mean, ymax = Mean + SD, ymin = Mean - SD, fill = pH ))
  k + geom_bar(stat="identity")  + geom_errorbar( width = 0.3) + scale_fill_brewer(type="qual", palette="Paired") + labs(title=paste(plot_title, add_title))
   ggsave(filename = paste(out_dir, plot_title, "_", add_title, "_genewise_FPKM", ".pdf", sep =""), width = 6, height = 6)
  }

# test plotting fpkm data
#plot_genewise_fpkm(gene_data, "titleA", "titleB")

#####
# function: plot_interaction
# makes an interaction plot for a given gene
#####

# plot_interaction <- function(gene_data, plot_title, add_title) {
#   h <- ggplot(gene_data, aes(x = pH, y = Mean, ymax = Mean + SD, ymin = Mean - SD, color=Mom))
#   h + geom_errorbar(color="black", width=0.1) + geom_point(aes(shape=Mom), color="black",size = 6) + geom_line(aes(group=Mom), size=0.5, linetype=1) + labs(title=paste(plot_title, add_title), x = "pH treatment", y = "mean FPKM") + scale_color_manual(values=c("#3399FF", "#66CC00", "#CC0066"))
#   ggsave(filename = paste(out_dir, plot_title, "_", add_title, "_interactionPlot_FPKM", ".pdf", sep =""), width = 6, height = 6)
# }

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
  plot_genewise_fpkm(gene_data, gene_name, graph_title)
  #plot_interaction(gene_data, gene_name, graph_title)
}

# generate genewise plots based on a gene of interest
GOI_data("comp100107_c0_seq1", transcripts, all_fpkm[,c(new.hatch.groups$Sample.ID)], "Carbonic Anhydrase")

###
#
#  Example: read in a list of genes of interest, and make plots
#
###

desired_graphs = read.csv(file.choose())  #a text file with two columns, the geneID (gene_name) and a descriptor (graph_title)

for (i in 1:nrow(desired_graphs)) {
  GOI_data(desired_graphs[i,1], transcripts, all_fpkm[,-c(1:8)], desired_graphs[i,2])
}

