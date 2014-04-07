
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


#####
#
# User-set paramaters
#
#####

outfile_prefix <- "outfile"

# set working directory
setwd("/Users/Shared/BigCB_Insect_Project/Dicosmo_project/Dicosmo_DGE")

# get the data files produced by import_express_data.R
load(file="diff_expression_data.RData")

# get list of transcripts that we want to cluster and make a heatmap from
targetList <- read.table("list_of_target_genes_for_heatmap.txt", header=TRUE)

# define number of clusters
k = 16

#####
#
# Generate clusters and heatmap
#
#####

data <- data.matrix(all_fpkm[(rownames(all_fpkm) %in% targetList$transcript_id),]) # make matrix of fpkm values using only those transcripts found in targetList

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
