#!/usr/bin/R

# Scott Fay 
# 2013 Oct 15
# updated 9 April 2014

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
setwd("/Users/Shared/RKC_project/RKC5/RKC5_new_mapping_ExpressOutput/RKC5_new_ExpressOutput")
# set output directory
out_dir <- "/Users/Shared/RKC_project/RKC5/RKC5_new_mapping_ExpressOutput/"
# get annotation file
# "trinotate_annotation_report.txt" is a Trinotate output excel file, exported as tab-delimited
transcripts <- read.csv("/Users/scottfay/Dropbox/RKC_ALL_Mapped_Annot/ALL_RKC_trinotate_annotation_report.csv", stringsAsFactors=FALSE)

# define groups
groups <- read.delim("/Users/Shared/RKC_project/RKC5/RKC5_new_mapping_ExpressOutput/RKC5_groups.txt", colClasses = c("character", rep("factor",4)))
str(groups)

##################
#
# reorder groups so that ambient pH comes first
#
#####################
groups2 = groups
groups2$pH = factor(groups$pH,levels(groups$pH)[c(3,2,1)])
groups2
str(groups2)
groups = groups2

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


# validate the sample identity:
#   print out the list of samples
#   ask the user to compare it to the sample/group definition
#   get permission to proceed
cat("\nPlease compare the list of directories to the list of samples/groups:\n\n")
data.frame(Directory=libs, groups)
cat("Does everything look copacetic?\nIf not edit groups file to match sample IDs before importing again.")

# Rename columns
colnames(all_count) <- groups$sample.ID
head(all_count)
colnames(all_fpkm) <- groups$sample.ID
head(all_fpkm)
# now the data is ready

#####
#
# filter out transcripts by mean count per gene
# NOTE: another kind of filtering, that takes into account library size, is used below.  Prefer that method over this one.
#
####
# 
# filterVal <- 10    # filters out mean count per gene of less than this value
# means <- rowMeans(all_count)
# filter <- means >= filterVal
#  table(filter)
#  filtered.all.count <- all_count[filter,]
#  dim(filtered.all.count)

#####
#
# Plot total number of mapped reads to known genes: look for outliers of library size
#
#####

##library(RColorBrewer)
#colors <- brewer.pal(9, "Set1")
boxplot(log2(all_count+1), las=2, main="After filter") # before filter
boxplot(log2(filtered.all.count+1), las=2, main=paste("After filtering: >", filterVal, " mean counts", sep="")) # after filter

#####
#
# Drop any poor quality libraries (sample 232 is funky on the boxplot (small))
#
#####

drops <- c("323")  # list the library names to be removed here
# 
# filtered.all.count <- filtered.all.count[,!(names(filtered.all.count) %in% drops)]
# 
all_fpkm <- all_fpkm[,!(names(all_fpkm) %in% drops)]
all_count <- all_count[,!(names(all_count) %in% drops)]

# head(filtered.all.count)
groups <- groups[!(groups$sample.ID %in% drops),]
# boxplot(log2(filtered.all.count+1), las=2, main=paste("After filtering: >", filterVal, " mean counts", sep="")) # visualize again

#####
#
# Make EdgeR objects and normalize
#
#####

x <- DGEList(counts=all_count)

# TMM normalization
x.tmm <- calcNormFactors(x, method="TMM")

# upper quartile normalization
x.uq <- calcNormFactors(x, method="upperquartile")

#####
#
# MDS plots
#
#####

plotMDS(x , main = "MDS Plot for Count Data, Raw", labels = colnames( x$counts ), cex=0.7)

plotMDS(x.tmm , main = "MDS Plot for Count Data, TMM Norm", labels = colnames( x$counts ), cex=0.7)

plotMDS(x.uq , main = "MDS Plot for Count Data, UpperQuartile Norm", labels = colnames( x$counts ), cex=0.7)

# other labels:
plotMDS(x.tmm , main = "MDS Plot for Count Data by Individual", labels = groups$mother, cex=0.7)
plotMDS(x.tmm , main = "MDS Plot for Count Data by Tissue", labels = groups$tissue, cex=0.7)
plotMDS(x.tmm , main = "MDS Plot for Count Data by pH", labels = groups$pH, cex=0.7)

# more dimensions, defined in dim.plot()
plotMDS(x.tmm , main = "MDS Plot for Count Data", labels = groups$tissue, , dim.plot=c(3,4), cex=0.7)

#####
#
# Partition data by Tissue: create four new data sets, one for each tissue
#
#####

for (tissue in levels(groups$tissue)) {
  assign(paste("groups.", tissue, sep=""), groups[groups$tissue == tissue,]) # makes a "groups" object for the particular tissue
  assign(paste("y.", tissue, sep=""), DGEList(all_count[,groups[groups$tissue == tissue, 1]])) # makes a DGE object for the particular tissue
}



#####
# 
# Test of differential expression
#
#####

# define the statistical model design matrix using the model.matrix function.  THis model is to look at effects of pH adjusting for the  differences among tissues. (from EdgeR documents design matrix in section 4.4.7)
str(groups)
pH <- groups$pH
mother <- groups$mother
tissue <- groups$tissue

design <- model.matrix(~tissue+pH) #the full data set, not broken down by tissues.  What we used.
design

#design <- model.matrix(~pH+pH:tissue) #the full data set, not broken down by tissues, partial interactive effects of pH and tissue
#design <- model.matrix(~pH*tissue) #the full data set, not broken down by tissues, partial interactive effects of pH and tissue

# FILTERING: filter out transcripts that do not have at least five reads out of 1,000,000 reads mapped in at least 4 samples
# NOTE: BCV values look much better when you do this kind of filtering than just mean > a certain threshhold.
x2.tmm <- x.tmm[rowSums(1e+06 * x.tmm$counts/expandAsMatrix(x.tmm$samples$lib.size, dim(x.tmm)) > 5) >= 4, ]

# reset the library sizes after filtering
x2.tmm$samples$lib.size <- colSums(x2.tmm$counts)

# estimate dispersion for GLM fit, common, trended, and tagwise
# see edgeR docs for references on dispersion estimates, e.g., "?estimateGLMCommonDisp"
x2.tmm <- estimateGLMCommonDisp(x2.tmm, design) 

# estimate trended dispersion for use in tagwise dispersion estimate
x2.tmm <- estimateGLMTrendedDisp(x2.tmm,design)
# estimate tagwise dispersion to be used in glmFit()
x2.tmm <- estimateGLMTagwiseDisp(x2.tmm,design)

# plot genewise biological coefficient of variation against gene abundance
plotBCV(x2.tmm)
# make plot as a pdf
pdf(file=paste(out_dir, "BCV_plot.pdf", sep=""), height=6, width=6)
plotBCV(x2.tmm, main = "Biological Coefficient of Variation")
dev.off()
# MDS plot again
plotMDS(x2.tmm , main = "MDS Plot for Count Data", labels = colnames( x2.tmm$counts ), cex=0.7)
# make a pdf
pdf(file=paste(out_dir, "MDS_plot.pdf", sep=""), height=6, width=6)
plotMDS(x2.tmm , main = "MDS Plot for Count Data", labels = colnames( y$counts ), cex=0.7)
dev.off()

# GLM
# fit the negative binomial generalized linear model (GLM) for each tag, creating a new fit object
fit <- glmFit(x2.tmm, design)
colnames(fit)

#####
# 
# Liklihood ratio test (LRT) to show the top genes
#
#####
#e.g., lrt6 is pH 7.5 vs. ambient (following example in EdgeR documentation section 3.3.2, pg 26)
lrt6 = glmLRT(fit, coef=6)
topTags(lrt6)
summary(decideTestsDGE(lrt6, p=0.05, adjust="BH"))
lrt5 = glmLRT(fit, coef=5)
topTags(lrt5)
summary(decideTestsDGE(lrt5, p=0.05, adjust="BH"))
lrt4 = glmLRT(fit, coef=4)
topTags(lrt4)
summary(decideTestsDGE(lrt4, p=0.05, adjust="BH"))
lrt3 = glmLRT(fit, coef=3)
topTags(lrt3)
summary(decideTestsDGE(lrt3, p=0.05, adjust="BH"))
lrt2 = glmLRT(fit, coef=2)
topTags(lrt2)
summary(decideTestsDGE(lrt2, p=0.05, adjust="BH"))

lrt7 <- glmLRT(fit, contrast=c(0,-1,1,0,0,0))
lrt8 <- glmLRT(fit, contrast=c(0,-1,0,1,0,0))
lrt9 <- glmLRT(fit, contrast=c(0,0,0,0,-1,1))
               
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

# call the get_DEGs() function for the relevant comparison made above, using the LRT fit as generated above...
DEGlist_lrt6 <- get_DEGs(lrt6, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt5 <- get_DEGs(lrt5, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt4 <- get_DEGs(lrt4, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt3 <- get_DEGs(lrt3, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt2 <- get_DEGs(lrt2, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt7 <- get_DEGs(lrt7, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt8 <- get_DEGs(lrt8, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt9 <- get_DEGs(lrt9, annot=transcripts, fdr=0.01, critFC=2)


#make lists with of the significantly differentially expressed putative transcripts

RKC5_DEGlist <- rbind(DEGlist_lrt2,DEGlist_lrt3,DEGlist_lrt4,DEGlist_lrt5,DEGlist_lrt6,DEGlist_lrt7,DEGlist_lrt8,DEGlist_lrt9)
DE_transcripts <- RKC5_DEGlist$transcript_id
length(DE_transcripts)
unique_DE_transcripts <- unique(DE_transcripts)
length(unique_DE_transcripts)  # number of unique differentially expressed transcripts across all comparisons
write.csv(unique_DE_transcripts, paste(out_dir, "all_comparison_DEG_list.csv", sep=""))

RKC5_DEGlist_only_tissue <- rbind(DEGlist_lrt2,DEGlist_lrt3,DEGlist_lrt4,DEGlist_lrt7,DEGlist_lrt8)
unique_tissue_DE_transcripts <- unique(RKC5_DEGlist_only_tissue$transcript_id)
length(unique_tissue_DE_transcripts) # number of uniquely differentially expressed transcripts across pH comparisons
write.csv(unique_tissue_DE_transcripts, paste(out_dir, "tissue_comparison_DEG_list.csv", sep=""))

RKC5_DEGlist_only_pH <- rbind(DEGlist_lrt5,DEGlist_lrt6,DEGlist_lrt9)
unique_pH_DE_transcripts <- unique(RKC5_DEGlist_only_pH$transcript_id)
length(unique_pH_DE_transcripts) # number of uniquely differentially expressed transcripts across pH comparisons
write.csv(unique_pH_DE_transcripts, paste(out_dir, "pH_comparison_DEG_list.csv", sep=""))



# read in the groupings again, this time with the columns as characters (not factors)
names <- read.delim("../RKC5_groups.txt", colClasses = c(rep("character",4)))
names = names[-9,]

# function to concatenate columns for names
concat_name <- function(x){
  paste(x[1],x[2],x[3],x[4],sep="_")
}

transcripts$CloneID <- transcripts$transcript_id
head(transcripts)


###
# all
###
# make matrix of fpkm values using only those transcripts found in all targetList
FPKM_all <- data.matrix(all_fpkm[(rownames(all_fpkm) %in% unique_DE_transcripts),])
# log2 transform, mean center rows
FPKM_all = log2(FPKM_all+1)
centered_data = t(scale(t(FPKM_all), scale=F)) # center rows, mean substracted
# make data frame of centered_data matrix
centeredDF <- as.data.frame(centered_data)
# concatenate column names from groupings to make the column names more informative
colnames(centeredDF) <- apply(names, 1, concat_name)
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
write.table(centeredDF, file=paste(out_dir, "centered_data_all_comparison_DEGs.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")


###
# tissue
###
# make matrix of fpkm values using only those transcripts found in tissue targetList
FPKM_tissue <- data.matrix(all_fpkm[(rownames(all_fpkm) %in% unique_tissue_DE_transcripts),])
# log2 transform, mean center rows
FPKM_tissue = log2(FPKM_tissue+1)
centered_data = t(scale(t(FPKM_tissue), scale=F)) # center rows, mean substracted
# make data frame of centered_data matrix
centeredDF <- as.data.frame(centered_data)
# concatenate column names from groupings to make the column names more informative
colnames(centeredDF) <- apply(names, 1, concat_name)
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
write.table(centeredDF, file=paste(out_dir, "centered_data_tissue_comparison_DEGs.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")


###
# pH
###
# make matrix of fpkm values using only those transcripts found in pH targetList
FPKM_pH <- data.matrix(all_fpkm[(rownames(all_fpkm) %in% unique_pH_DE_transcripts),])
# log2 transform, mean center rows
FPKM_pH = log2(FPKM_pH+1)
centered_data = t(scale(t(FPKM_pH), scale=F)) # center rows, mean substracted
# make data frame of centered_data matrix
centeredDF <- as.data.frame(centered_data)
# concatenate column names from groupings to make the column names more informative
colnames(centeredDF) <- apply(names, 1, concat_name)
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
write.table(centeredDF, file=paste(out_dir, "centered_data_pH_comparison_DEGs.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")



##
#
# Get all the transcripts that were being used (for annotation of less than all the data)
#
##
head(FPKM_all)
FPKM_all_Table=as.table(FPKM_all)
ALL_RKC_transcripts_RKC5 = rownames(FPKM_all_Table)
head(ALL_RKC_transcripts_RKC5)
ALL_RKC_transcripts_RKC5 = as.data.frame(ALL_RKC_transcripts_RKC5) #700 obs
write.table(ALL_RKC_transcripts_RKC5,"RKC5_All_RKC_transcriptID.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

###
#
# Check to see if ALL_RKC_transcripts_RKC5 contains transcripts not already included in ALL_RKC_transcripts from RKC2.1
#
##

RKC2.1_transcript = read.table(file.choose())  # read in the file containing the transcript IDs identified in RKC2.1 (*32955 data)
head(RKC2.1_transcript)

# see if any transcript IDs from ALL_RKC_transcripts_RKC5 are not in RKC2.1_transcript
colnames(ALL_RKC_transcripts_RKC5)=c("V1")
RKC5_2.1_Transcripts <- rbind(RKC2.1_transcript,ALL_RKC_transcripts_RKC5)
unique_RKC5_2.1_Transcripts <- unique(RKC5_2.1_Transcripts) #33243 data (means that there are 288 transcripts from RKC5 not in RKC2.1)
# what are those transcripts?  use grep in the bash shell to find out.
# grep -v 

testXX = for (i in 1:700) {
  grep(ALL_RKC_transcripts_RKC5[i,], RKC5_2.1_Transcripts, invert=TRUE)
  #tempframe <- read.delim(paste('./', dirname, '/results.xprs', sep=''))
  #assign(dirname, tempframe[order(tempframe$target_id),]) # generate dataframe, organizing rows by target_id
}


#####
#
# Normalize EdgeR objects by tissue  (NOTE NOT working properly...?? why?).  Didn't do this yet.
#
#####

# TMM normalization
y.antG.tmm <- calcNormFactors(y.antG, method="TMM")
y.postG.tmm <- calcNormFactors(y.postG, method="TMM")
y.cuti.tmm <- calcNormFactors(y.cuti, method="TMM")
y.heart.tmm <- calcNormFactors(y.heart, method="TMM")

# upper quartile normalization
y.antG.uq <- calcNormFactors(y.antG, method="upperquartile")
y.postG.uq <- calcNormFactors(y.postG, method="upperquartile")
y.cuti.uq <- calcNormFactors(y.cuti, method="upperquartile")
y.heart.uq <- calcNormFactors(y.heart, method="upperquartile")


#######
#
# Analyze just pH effects in antG
#
#######

#####
# 
# Test of differential expression
#
#####

# define the statistical model design matrix using the model.matrix function.  THis model is to look at effects of pH adjusting for the  differences among tissues. (from EdgeR documents design matrix in section 4.4.7)
str(groups)
pH <- groups$pH
mother <- groups$mother

design <- model.matrix(~pH)
design

#design <- model.matrix(~pH+pH:tissue) #the full data set, not broken down by tissues, partial interactive effects of pH and tissue
#design <- model.matrix(~pH*tissue) #the full data set, not broken down by tissues, partial interactive effects of pH and tissue

# FILTERING: filter out transcripts that do not have at least five reads out of 1,000,000 reads mapped in at least 4 samples
# NOTE: BCV values look much better when you do this kind of filtering than just mean > a certain threshhold.
y2antG.tmm <- y.antG[rowSums(1e+06 * y.antG$counts/expandAsMatrix(y.antG$samples$lib.size, dim(y.antG)) > 5) >= 4, ]

# reset the library sizes after filtering
y2antG.tmm$samples$lib.size <- colSums(y2antG.tmm$counts)

# estimate dispersion for GLM fit, common, trended, and tagwise
# see edgeR docs for references on dispersion estimates, e.g., "?estimateGLMCommonDisp"
y2antG.tmm <- estimateGLMCommonDisp(y2antG.tmm, design) 

# estimate trended dispersion for use in tagwise dispersion estimate
y2.antG.tmm <- estimateGLMTrendedDisp(y2.antG.tmm,design)
# estimate tagwise dispersion to be used in glmFit()
y2.antG.tmm <- estimateGLMTagwiseDisp(y2.antG.tmm,design)

# plot genewise biological coefficient of variation against gene abundance
plotBCV(x2.tmm)
# make plot as a pdf
pdf(file=paste(out_dir, "BCV_plot.pdf", sep=""), height=6, width=6)
plotBCV(x2.tmm, main = "Biological Coefficient of Variation")
dev.off()
# MDS plot again
plotMDS(x2.tmm , main = "MDS Plot for Count Data", labels = colnames( x2.tmm$counts ), cex=0.7)
# make a pdf
pdf(file=paste(out_dir, "MDS_plot.pdf", sep=""), height=6, width=6)
plotMDS(x2.tmm , main = "MDS Plot for Count Data", labels = colnames( y$counts ), cex=0.7)
dev.off()

# GLM
# fit the negative binomial generalized linear model (GLM) for each tag, creating a new fit object
fit <- glmFit(x2.tmm, design)
colnames(fit)

#####
# 
# Liklihood ratio test (LRT) to show the top genes
#
#####
#e.g., lrt6 is pH 7.5 vs. ambient (following example in EdgeR documentation section 3.3.2, pg 26)
lrt6 = glmLRT(fit, coef=6)
topTags(lrt6)
summary(decideTestsDGE(lrt6, p=0.05, adjust="BH"))
lrt5 = glmLRT(fit, coef=5)
topTags(lrt5)
summary(decideTestsDGE(lrt5, p=0.05, adjust="BH"))
lrt4 = glmLRT(fit, coef=4)
topTags(lrt4)
summary(decideTestsDGE(lrt4, p=0.05, adjust="BH"))
lrt3 = glmLRT(fit, coef=3)
topTags(lrt3)
summary(decideTestsDGE(lrt3, p=0.05, adjust="BH"))
lrt2 = glmLRT(fit, coef=2)
topTags(lrt2)
summary(decideTestsDGE(lrt2, p=0.05, adjust="BH"))

lrt7 <- glmLRT(fit, contrast=c(0,-1,1,0,0,0))
lrt8 <- glmLRT(fit, contrast=c(0,-1,0,1,0,0))
lrt9 <- glmLRT(fit, contrast=c(0,0,0,0,-1,1))

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

# call the get_DEGs() function for the relevant comparison made above, using the LRT fit as generated above...
DEGlist_lrt6 <- get_DEGs(lrt6, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt5 <- get_DEGs(lrt5, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt4 <- get_DEGs(lrt4, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt3 <- get_DEGs(lrt3, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt2 <- get_DEGs(lrt2, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt7 <- get_DEGs(lrt7, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt8 <- get_DEGs(lrt8, annot=transcripts, fdr=0.01, critFC=2)
DEGlist_lrt9 <- get_DEGs(lrt9, annot=transcripts, fdr=0.01, critFC=2)


#make lists with of the significantly differentially expressed putative transcripts

RKC5_DEGlist <- rbind(DEGlist_lrt2,DEGlist_lrt3,DEGlist_lrt4,DEGlist_lrt5,DEGlist_lrt6,DEGlist_lrt7,DEGlist_lrt8,DEGlist_lrt9)
DE_transcripts <- RKC5_DEGlist$transcript_id
length(DE_transcripts)
unique_DE_transcripts <- unique(DE_transcripts)
length(unique_DE_transcripts)  # number of unique differentially expressed transcripts across all comparisons
write.csv(unique_DE_transcripts, paste(out_dir, "all_comparison_DEG_list.csv", sep=""))

RKC5_DEGlist_only_tissue <- rbind(DEGlist_lrt2,DEGlist_lrt3,DEGlist_lrt4,DEGlist_lrt7,DEGlist_lrt8)
unique_tissue_DE_transcripts <- unique(RKC5_DEGlist_only_tissue$transcript_id)
length(unique_tissue_DE_transcripts) # number of uniquely differentially expressed transcripts across pH comparisons
write.csv(unique_tissue_DE_transcripts, paste(out_dir, "tissue_comparison_DEG_list.csv", sep=""))

RKC5_DEGlist_only_pH <- rbind(DEGlist_lrt5,DEGlist_lrt6,DEGlist_lrt9)
unique_pH_DE_transcripts <- unique(RKC5_DEGlist_only_pH$transcript_id)
length(unique_pH_DE_transcripts) # number of uniquely differentially expressed transcripts across pH comparisons
write.csv(unique_pH_DE_transcripts, paste(out_dir, "pH_comparison_DEG_list.csv", sep=""))



# read in the groupings again, this time with the columns as characters (not factors)
names <- read.delim("../RKC5_groups.txt", colClasses = c(rep("character",4)))
names = names[-9,]

# function to concatenate columns for names
concat_name <- function(x){
  paste(x[1],x[2],x[3],x[4],sep="_")
}

transcripts$CloneID <- transcripts$transcript_id
head(transcripts)


###
# all
###
# make matrix of fpkm values using only those transcripts found in all targetList
FPKM_all <- data.matrix(all_fpkm[(rownames(all_fpkm) %in% unique_DE_transcripts),])
# log2 transform, mean center rows
FPKM_all = log2(FPKM_all+1)
centered_data = t(scale(t(FPKM_all), scale=F)) # center rows, mean substracted
# make data frame of centered_data matrix
centeredDF <- as.data.frame(centered_data)
# concatenate column names from groupings to make the column names more informative
colnames(centeredDF) <- apply(names, 1, concat_name)
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
write.table(centeredDF, file=paste(out_dir, "centered_data_all_comparison_DEGs.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")

