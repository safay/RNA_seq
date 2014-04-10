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



#####
#
# JS plyr version
#
#####
pH # use as.factor from group file, or make this as a factor
temp # use as.factor from group file, or make this as a factor

library(plyr)
library(ggplot2)

gene="comp0_c0_seq1"

get_gene_data <- function(gene, all_fpkm) {
  gene_data <- data.frame(fpkm=as.numeric(all_fpkm[gene,]), pH=pH, Temperature=temp)
  fpkm_mean<-ddply(gene_data,.(pH,Temperature), function(d) mean(d$fpkm)) # if group file has columns called pH and Temperature
  SD<-ddply(gene_data,.(pH,Temperature), function(d) sd(d$fpkm)) # if group file has columns called pH and Temperature
  names(fpkm_mean)=c("pH","Temperature","Mean") # rename
  names(SD)=c("pH","Temperature","SD") # rename
  out_frame=merge(fpkm_mean,SD) #make one data frame
  out_frame=out_frame[,c(3,4,1,2)] # reorder columns
  return(out_frame)
}

gene_data <- get_gene_data("comp0_c0_seq1", all_fpkm)  #the gene in all_fpkm for which data will be returned





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
  h <- ggplot(gene_data, aes(x = pH, y = Mean, ymax = Mean + SD, ymin = Mean - SD, color=Temperature))
  h + geom_errorbar(color="black", width=0.1) + geom_point(aes(shape=Temperature), color="black",size = 6) + geom_line(aes(group=Temperature), size=0.5, linetype=1) + labs(title=paste(plot_title, add_title), x = "pH treatment", y = "mean FPKM") + scale_color_manual(values=c("#3399FF", "#66CC00", "#CC0066"))
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

for (i in 1:nrow(desired_graphs)) {
  GOI_data(desired_graphs[i,1], transcripts, all_fpkm, desired_graphs[i,2])
}



