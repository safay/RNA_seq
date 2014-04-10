#!/usr/bin/Rscript

#####
#
# Using data from "import_express_data.R,"
# prepares input file for Cluster 3.0, open-source clustering software
#   http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm
# 
# Usage:
#   command line:
#   prep_clustering.R 
#
#####

#choose the output "diff_expression_data.RData" from get_express_data.R here to get the file path
getfilepath <- file.choose()

# set the working directory with the output from get_express_data.R
setwd(dirname(getfilepath))

# load FPKM values, counts and groupings (output of import_express_data.R script)
load(getfilepath)
group

# get list of transcripts that we want to cluster and make a heatmap from
targetList <- read.csv(file.choose())

#####
#
# Transform
#
#####

# make matrix of fpkm values using only those transcripts found in targetList
data <- data.matrix(all_fpkm[(rownames(all_fpkm) %in% targetList$transcript_id),])

# log2 transform, mean center rows
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) # center rows, mean substracted

### Save centered_data
write.csv(centered_data, "centered_data.csv")

# this file ready to cluster with Cluster 3.0

DEGlist_lrt_inter_0

# read in the groupings again, this time with the columns as characters (not factors)
names <- read.delim("../RKC5_groups.txt", colClasses = c(rep("character",4)))
names = names[-9,]

# function to concatenate columns for names
concat_name <- function(x){
  paste(x[1],x[2],x[3],x[4],sep="_")
}

# concatenate column names from groupings to make the column names more informative
new_column_names <- apply(names, 1, concat_name)

transcripts$CloneID <- transcripts$transcript_id
head(transcripts)

#####
#
# function to write cluster file
#   arguments:
#     all_fpkm = all_fpkm, a dataframe of fpkm values
#     unique_DE_transcripts, a list of differentially expressed genes
#     outfile = path and name of output file
#     cnames = vector of new column names
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
  write.table(centeredDF, file=paste(out_dir, "centered_data_all_comparison_DEGs.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")
}
