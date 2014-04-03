#!/usr/bin/R

# Scott Fay 
# 15 Jan 2013

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

# Load in FPKM object made in get_express_data.R

# load group file; should be in the same directory as FPKM object 
group <- read.delim(" PATH TO GROUP FILE; FILENAME")

# get annotation file
# define this on the command line
# "trinotate_annotation_report.txt" is a Trinotate output excel file, exported as tab-delimited
transcripts <- read.delim("/Users/scottfay/annotation/Calineuria/trinotate_annotation_report.txt")

# get the transcript ID from the command line too

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
  mean_15 <- apply(fpkm[gene,1:5], 1, mean)
  mean_20 <- apply(fpkm[gene,6:10], 1, mean)
  mean_25 <- apply(fpkm[gene,11:15], 1, mean)
  mean_30 <- apply(fpkm[gene,16:20], 1, mean)
  sd_15 <- apply(fpkm[gene,1:5], 1, sd)
  sd_20 <- apply(fpkm[gene,6:10], 1, sd)
  sd_25 <- apply(fpkm[gene,11:15], 1, sd)
  sd_30 <- apply(fpkm[gene,16:20], 1, sd)
  fpkm_mean <- c(mean_15, mean_20, mean_25, mean_30)
  SD <- c(sd_15, sd_20, sd_25, sd_30)
  temp <- c("15C", "20C", "25C", "30C")
  out_frame <- data.frame(fpkm_mean, SD, temp)
  return(out_frame)
}

# test getting fpkm data for a specific gene
#gene_data <- get_gene_data("comp63351_c0_seq3", all_fpkm)

#####
# Function: get_annot
# gets annotation for a particular gene
#####

get_annot <- function(gene_in, annotation) {
  L <- annotation$transcript_id == gene_in
  cat("transcript: ", as.vector(annotation[L,]$transcript_id), "\n")
  cat("top blast hit: ", as.vector(annotation[L,]$Top_BLASTX_hit), "\n")
  cat("gene ontology: ", as.vector(annotation[L,]$gene_ontology), "\n")
  cat("PFam: ", as.vector(annotation[L,]$Pfam), "\n")
  cat("Signal Peptide: ", as.vector(annotation[L,]$SignalP), "\n")
  cat("transmembrane domain: ", as.vector(annotation[L,]$TmHMM), "\n")
  cat("eggnog: ", as.vector(annotation[L,]$eggnog), "\n")
}

# test getting annotation data
#get_annot("comp61673_c0_seq1", transcripts)

#####
# Function: plot_genewise_fpkm
# plots fpkm for the specified gene
#####

plot_genewise_fpkm <- function(gene_in, plot_title, add_title) {
  k <- ggplot(gene_in, aes( x = temp, y = fpkm_mean, ymax = fpkm_mean + SD, ymin = fpkm_mean - SD ) )
  k + geom_point( size = 3 ) + geom_errorbar(width = 0.4) + theme_bw() + theme(plot.title = element_text(size = rel(0.9)), axis.title = element_text(size = rel(0.8))) + labs(title=paste(plot_title, "\n", add_title), x = "Temperature", y = "mean FPKM")
  ggsave(filename = paste(out_dir, plot_title, "_", add_title, "_genewise_FPKM", ".pdf", sep =""), width = 2, height = 2.5)
}

#plot_genewise_fpkm(gene_data, "comp63351_c0_seq3", " expression")

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
}

# generate genewise plots based on a gene of interest
GOI_data("comp55133_c0_seq1", transcripts, all_fpkm, "Zinc Finger")

