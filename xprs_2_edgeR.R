#!/usr/bin/R

# Scott Fay 
# 2013 Oct 15

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
setwd("/Users/Shared/BigCB_Insect_Project/Ptero_DGE/Pt_DGE_noSpike/Ptero_DGE")
# set output directory
out_dir <- "/Users/Shared/BigCB_Insect_Project/Ptero_DGE/Pt_DGE_noSpike/DGE_output"
# get annotation file
# "trinotate_annotation_report.txt" is a Trinotate output excel file, exported as tab-delimited
transcripts <- read.delim("/Users/Shared/BigCB_Insect_Project/Ptero_DGE/trinotate_annotation_report.txt")
# define groups
group <- factor(c(rep("15", 5), rep("20",5), rep("25",5), rep("30",5)))

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

# now the data is ready

# remove spike option
# retreive the spike IDs generated using blast
#spikes <- scan("~/Google\ Drive/RNA_seq/Data/RKC_Assemblies/spike_transcript_IDs", what = "character")
# remove the spikes transcripts from the tables
#all_count <- all_count[!rownames(all_count) %in% spikes, ]
#all_fpkm <- all_fpkm[!rownames(all_fpkm) %in% spikes, ]

##########
#
# Set up differential gene expression (DGE) analysis, generate summary plots:
#   define a model matrix *
#   estimate dispersion
#   make a plot of biological coefficient of variation vs. log(CPM)
#   make a multidimensional scaling plot of total counts
#   fit GLM for each transcript
#     * Note, how to define the statistical model is completely dependent on your experimental design and what question(s) you want to address.  See the EdgeR documentation as a first-pass guide: http://www.bioconductor.org/packages/2.12/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#
##########


# DGEList makes an EdgeR object
y <- DGEList(counts=all_count, group=group)

# define the statistical model, a design matrix using the model.matrix function
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
colnames(design)
design
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
#   lrt = glmLRT object of the comparison you're making
#   annot = transcripts annotation data frame
#   fdr = false discovery rate (default 0.05)
#   critFC, a threshold fold change (default 0, i.e., all)
#####

get_DEGs <- function(lrt, annot, fdr=0.05, critFC=0, onlyAnnot=FALSE) {
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
  cat("percent all transcripts with Pfam or BlastX or BlastP annotation:", sum((annot$Pfam != ".") | (annot$Top_BLASTX_hit != ".") | (annot$Top_BLASTP_hit != ".")) * 100 / nrow(annot), "\n")
  cat("percent toptags with Pfam or BlastX or BlastP annotation:", sum((join$Pfam != ".") | (join$Top_BLASTP_hit != ".") | (join$Top_BLASTX_hit != ".")) * 100 / nrow(join), "\n")
  cat("number of genes upregulated, (FC > 2):", sum(join$logFC > 1), "\n")
  cat("number of genes downregulated, (FC < 0.5):", sum(join$logFC < -1 ), "\n")
  cat("...saving .csv file with all topTags:", paste( out_dir, "topTags_", tTags$comparison, "_", fdr, "_FDR.csv", sep="" ), "\n" )
  write.csv( join, file=paste(out_dir, "topTags_", tTags$comparison, "_", fdr, "_FDR.csv", sep="" ) )
  cat("...saving .csv file with pFam OR TopBlastX/P-Hit annotated topTags:", paste(out_dir, "topTags_", tTags$comparison, "_", fdr, "_FDR_only_w_annot.csv", sep=""), "\n" )
  # filter based on the threshold fold change, critFC
  join <- join[ abs(join$logFC) > sqrt(critFC) , ]
  # write a csv and return the relevant dataframe
  if(onlyAnnot) {
    cat("saving with only annot")
    write.csv( join[ (join$Pfam != ".") | (join$Top_BLASTX_hit != ".") | (join$Top_BLASTP_hit != ".") , ], file=paste(out_dir, "topTags_", tTags$comparison, "_", fdr, "_FDR_only_w_Annot.csv", sep="" ) )
    return(join[ (join$Pfam != ".") | (join$Top_BLASTX_hit != ".") | (join$Top_BLASTP_hit != ".") , ])
  } else {
    cat("saving toptags")
    write.csv( join, file=paste(out_dir, "topTags_", tTags$comparison, "_", fdr, "_FDR.csv", sep="" ) )
    return(join)    
  }
}

# perform the desired LRT based on what comparison/contrast you want to make
# Pt_15vs30 <- glmLRT(fit, coef=1)
Pt_15vs30 <- glmLRT(fit, contrast=c(-1,0,0,1))

# call the get_DEGs() function for the relevant comparison made above, using the LRT fit as generated above...
xxDEGlist_Pt_15vs30_annot_FC_10 <- get_DEGs(Pt_15vs30, annot=transcripts, fdr=0.00001, critFC=10, onlyAnnot=TRUE)

# sort by logFC
#sorted_tags <- DEGlist_Pt_15vs30[with(DEGlist_Pt_15vs30, order(logFC)), ]

# write a .csv based on a threshold FC, in this case, ( FC < 0.25 ) and ( FC > 4 )
write.csv(DEGlist_Pt_15vs30[ abs(DEGlist_Pt_15vs30$logFC) > 2 , ], file=paste( out_dir, "topTags_DEGlist_15vs30_FDR0_00001_PfamHits_logFC_greater_than_2.csv", sep="" ))

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
#####

GOI_data <- function(gene_name, annotation, fpkm_table, graph_title) {
  gene_data <- get_gene_data(gene_name, fpkm_table)
  get_annot(gene_name, annotation)
  plot_genewise_fpkm(gene_data, gene_name, graph_title)
  plot_interaction(gene_data, gene_name, graph_title)
}

# generate genewise plots based on a gene of interest
GOI_data("comp100107_c0_seq1", transcripts, all_fpkm, "Carbonic Anhydrase")


#####
#
# Clustering and heatmap
#
#####

# Adapted from Trinity differential gene expression analysis found here:
# http://trinityrnaseq.sourceforge.net/analysis/diff_expression_analysis.html
# from analyze_diff_expr.pl

library(cluster)
library(gplots)
library(Biobase)

outfile_prefix <- "outfile"

tTagList <- FC_10_Pt_15vs30

# define number of clusters
k = 16

data <- data.matrix(all_fpkm[(rownames(all_fpkm) %in% tTagList$transcript_id),]) # make matrix of fpkm values using only those transcripts found in tTagList

## generate correlation matrix
cr = cor(data, method='spearman')
## log2 transform, mean center rows
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) # center rows, mean substracted
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
myheatcol = redgreen(75)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=k)
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
#postscript(file=paste(out_dir,"diff_expr_matrix_file.heatmap.eps",sep=""), horizontal=FALSE, width=18, height=8, paper="special")
heatmap.2(centered_data, dendrogram="both", Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=2.5, margins=c(15,15), lhei=c(0.3,2), lwid=c(2.5,4))
#dev.off()

# prep for clustering
max_cluster_count = max(gene_partition_assignments)
gene_names = rownames(data)
num_cols = length(data[1,])

# partition data into subclusters, print plots
for (i in 1:max_cluster_count) {
  partition_i = (gene_partition_assignments == i)
  partition_data = centered_data[partition_i,]
  # if the partition involves only one row, then it returns a vector instead of a table\n";
  if (sum(partition_i) == 1) {
    dim(partition_data) = c(1,num_cols)
    colnames(partition_data) = colnames(centered_data)
    rownames(partition_data) = gene_names[partition_i]
  }
  assign(paste("subcluster_", i, sep=""), as.data.frame(partition_data))
  cluster_data <- as.data.frame(partition_data)
  cluster_fpkm_means <- rbind(data.frame(temp='15C', mean_fpkm=apply(cluster_data[,1:5], 1, mean), gene_id=rownames(cluster_data)), data.frame(temp='20C', mean_fpkm=apply(cluster_data[,6:10], 1, mean), gene_id=rownames(cluster_data)), data.frame(temp='25C', mean_fpkm=apply(cluster_data[,11:15], 1, mean), gene_id=rownames(cluster_data)), data.frame(temp='30C', mean_fpkm=apply(cluster_data[,16:20], 1, mean), gene_id=rownames(cluster_data)))
#  p <- ggplot(cluster_fpkm_means, aes(x=temp, y=mean_fpkm, group=gene_id))
  print(ggplot(cluster_fpkm_means, aes(x=temp, y=mean_fpkm, group=gene_id)) + geom_line(color=partition_colors[i],size=1, alpha=0.2) + geom_point(size = 3, alpha = 0.2) + theme_bw(base_size = 24) + ggtitle(NULL)) # + ggtitle(paste("Cluster_", i, sep="")) + ylab("median-centered log2(FPKM+1)") ))
}
