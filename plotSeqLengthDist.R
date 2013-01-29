#!/usr/bin/env Rscript 

# This R script, "plotSeqLengthDist," plots distributions of sequence lengths from an arbitrary number of FASTA files.
#
# command line usage:
#
# $ ./plotSeqLengthDist.R <file1 file2 ... > [-o <outputFile.ext>] [-s, -l, -d]
#
# -o : define an output file name; the extension given here defines the resulting form of output: pdf, jpeg, png
#
# use any combination of these three switches; each will create a separate plot
# -s : make stacked histogram
# -l : make superimposed (layered) histogram
# -d : make density plot
#
# Packages needed: ggplot2 and Biostrings (via Bioconductor)
# *** In the future, automated check for dependencies?
#
# To install ggplot2 in your R interpreter:
# > install.packages(“ggplot2”)
#
# To install Biostrings in your R interpreter:
# > source("http://bioconductor.org/biocLite.R")
# > biocLite("Biostrings")
#
# plotSeqLengthDist is licensed under Creative Commons, Attribution+Noncommercial+ShareAlike
# http://creativecommons.org/licenses/by-nc-sa/2.5/
# Author: Scott Fay
# scott.a.fay <at> gmail dot com
# 5 December 2012



# load required libraries
library(Biostrings)
library(ggplot2)

# define some functions...
# load individual files and, using biostrings, return a dataframe of two columns: transcript name and sequence lengths
loadFrame <- function(filepath, fileName)
{
  cat('loading ', fileName, '\n')
  seqs <- readDNAStringSet(filepath)               # Make a DNAStringSet object
  seqs_frame <- data.frame(length=width(seqs))     # Make a data frame, filling one column with sequence lengths
  seqs_frame["name"] <- fileName                   # Fill next column with the name
  return(seqs_frame)
}

# build up a combined data frame from each input file
makeCombinedFrame <- function() {
  allFrame <- data.frame()
  i <- 1                                    # initialize an index
  while (substr(args[i],1,1) != '-') {      # while command line arg doesn't have a hyphen at the start...
    if (file_test('-f',args[i])) {          # check if it's a file (***can this be changed to check if it is a FASTA file?)
      allFrame <- rbind(allFrame, loadFrame(args[i], basename(args[i])))         # load current frame into one combined frame
    }
    else {
      cat('error, check input file names and path \n')      # break if it's not a valid file
      break
    }
    if (is.na(args[i+1])) {
      cat('using default output filename, seqLengthDistOutput.pdf \n')
      break
    }
    i <- i+1                              # increment the index
  }
  return(allFrame)
}

# Make the plots using ggplot2
stackedH <- function(inputFrame) {
  cat('plotting stacked histogram \n')
  p <- ggplot(inputFrame, aes(x=length, fill=name)) + scale_x_log10() + geom_histogram(alpha=0.4, binwidth=0.01) + theme_bw(base_size = 18, base_family = "") + annotation_logticks(sides = 'b') + xlab('length (bp)') + guides(fill=guide_legend(title=NULL))
  ggsave(p, file=(paste(dirname(outFile), '/stacked_', basename(outFile), sep = '')), width=12, height=7)
}

layeredH <- function(inputFrame) {
  cat('plotting layered histogram \n')
  m <- unique(inputFrame$name)
  q <- ggplot(inputFrame, aes(x=length, fill=name)) + scale_x_log10() + theme_bw(base_size = 18, base_family = "") + annotation_logticks(sides = 'b') + xlab('length (bp)') + guides(fill=guide_legend(title=NULL))
  for (x in m) {
    q <- q + geom_histogram(data=subset(inputFrame,name==x), alpha=0.2, binwidth=0.01)
  }
  ggsave(q, file=(paste(dirname(outFile), '/layered_', basename(outFile), sep = '')), width=12, height=7)
}

densityP <- function(inputFrame) {
  cat('making density plot \n')
  r <- ggplot(inputFrame, aes(x=length, fill=name)) + geom_density(alpha=0.2) + scale_x_log10() + theme_bw(base_size = 18, base_family = "") + annotation_logticks(sides = 'b') + xlab('length (bp)') + guides(fill=guide_legend(title=NULL))
  ggsave(r, file=(paste(dirname(outFile), '/density_', basename(outFile), sep = '')), width=12, height=7)
}

# parse command line args to call each of the distribution plotting functions.
plotMaster <- function() {
  for (j in args) {
    if (j=='-s') stackedH(comboFrame) 
    else if (j=='-l') layeredH(comboFrame)
    else if (j=='-d') densityP(comboFrame)
  }
}

# Okay, get down to business
# Read in arguments from the command line
args=(commandArgs(TRUE))                

# Initialize default output file path and name; redefine it if the user wants
outFile <- './plotSeqLengthDistOutput.png'  
for (b in args) {
  if (b=='-o') {
    outFile <- args[match(b, args)+1]
  }
}

# Get the data
comboFrame <- makeCombinedFrame()

# Make the plots
plotMaster()
