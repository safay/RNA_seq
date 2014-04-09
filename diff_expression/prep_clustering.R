#!/usr/bin/Rscript

library(optparse)

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





