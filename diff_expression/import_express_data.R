#!/usr/bin/Rscript

# Scott Fay 
# 19 Mar 2014

#####
#
#
# Takes two files and prepares the data for futher differential gene expression analysis:
#  - output from eXpress (http://bio.math.berkeley.edu/eXpress/index.html)
#  - a tab-delimited table defining sample names and factors
#
# EXAMPLE USAGE:
# get_express_data.R </path/to/express/output/directories/> </path/groups_filename>
#  - First argument is a directory containing each of the xprs.output directories.
#  - Second argument is the path and filename to a "groups" file.
# Saves output to the directory containing the "groups" file.
#
# An example groups file, tab-delimited, with four factors, "mother", "pH", and "tissue."  
# Sample IDs are in the order of the directories containing the express output.
#
# sample.ID  mother	pH	tissue
# 309	1566	7.5	heart
# 313	1566	7.5	postG
# 314	1566	7.5	antG
# 315	1566	7.5	cuti
# 316	1575	7.5	heart
# 320	1575	7.5	postG
# 321	1575	7.5	antG
# 322	1575	7.5	cuti
# 323	1593	7.5	heart
# 327	1593	7.5	postG
# 328	1593	7.5	antG
# 329	1593	7.5	cuti
# 330	1609	7.8	heart
# 334	1609	7.8	postG
# 335	1609	7.8	antG
# 336	1609	7.8	cuti
# 337	1590	7.8	heart
# 341	1590	7.8	postG
# 342	1590	7.8	antG
# 343	1590	7.8	cuti
# 344	1592	7.8	heart
# 348	1592	7.8	postG
# 349	1592	7.8	antG
# 350	1592	7.8	cuti
# 365	1603	ambient	heart
# 369	1603	ambient	postG
# 370	1603	ambient	antG
# 371	1603	ambient	cuti
# 372	1572	ambient	heart
# 376	1572	ambient	postG
# 377	1572	ambient	antG
# 378	1572	ambient	cuti
# 379	1594	ambient	heart
# 383	1594	ambient	postG
# 384	1594	ambient	antG
# 385	1594	ambient	cuti
#
#####


# get filenames and paths from the command line
args = commandArgs(TRUE)
xprsDir <- args[1]
groupFile <- args[2]

# set working directory; must only contain directories containing eXpress output, no other files or directories
setwd(xprsDir)

# sets the groupfile directory as the output directory
out_dir <- dirname(groupFile)

# get a file defining groups
group <- read.delim(groupFile)
group

# creates a vector of the sequencing library directory names, "libs"
libs <- list.files()

# validate the sample identity:
#   print out the list of samples
#   ask the user to compare it to the sample/group definition
#   get permission to proceed
cat("\nPlease compare the list of directories to the list of samples/groups:\n\n")
data.frame(Directory=libs, group)
cat("\nDoes everything look copacetic?\ny / n ?\n")
response <- readLines(con="stdin",n=1)
while(TRUE) {  
  if (response == "n") {
    cat("Please edit the groups file to have the same order as the directories and try again...\n")
    quit(save="no")
    } else if (response == "y") {
      break
    } else {
      cat("y / n ?\n")
      response <- readLines(con="stdin",n=1)
    }
}

##########
#
# Import eXpress expression estimates as counts and FPKM values
#
##########

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

# clean up the dataframes a bit:
# create row names from target_id
rownames(all_count) <- all_count$target_id
rownames(all_fpkm) <- all_fpkm$target_id
# remove target_id column
all_count$target_id <- NULL
all_fpkm$target_id <- NULL
# change the column names to be the sample names from the grouping file
colnames(all_count) <- group$sample
colnames(all_fpkm) <- group$sample

#save the two new dataframes
save(all_count, all_fpkm, group, file = paste(out_dir, "/diff_expression_data.RData", sep=""))
