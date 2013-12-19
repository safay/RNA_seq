# this script gets data from the DGE .R script and Trinotate output ready for gene set enrichment analysis using the topGO R package
# http://www.bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
# basically this replaces the readMappings() function to take in Trinotate data
# see section 4.3 of the topGO documentation, "custom annotations"

# to install packages:
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("topGO")
# biocLite("Rgraphviz")

#####
# user-defined variables
#####

# a description of the experiment
expDescription <- "Calineuria thermal stress experiment"
# working directory with topTags .CSVs created by xprs_2_edgeR.R DGE script
setwd("/Users/Shared/BigCB_Insect_Project/Cali_project/Cali_DGE/")
# creates a vector of those DGE output files
tTag_output <- list.files()
tTag_output <- tTag_output[grepl("topTags.*csv", tTag_output)] # handy little grepl() function!
# annotation filename and path, tab-delimited file of Trinotate output, generated in MS Excel
annotFile <- "/Users/scottfay/annotation/Calineuria/trinotate_annotation_report.txt"
# number of enriched Gene Ontogeny terms to report
numNodes = 40
# where to put the output file
out_dir <- "/Users/Shared/BigCB_Insect_Project/Cali_project/Cali_DGE/"
# name of the output file
out_file <- "GSEA_Calineuria.txt"


#####
# Make full GO list from Trinotate-annotated transcriptome
#####
# get annotation file
annot <- read.delim(annotFile)
# create a list of GO terms
GO_list <- as.list(as.character(annot[['transcript_id']])) # generates list; element names are transcript IDs
GO_list <- as.list(setNames(as.character(annot$gene_ontology), as.character(annot$transcript_id))) # adds Gene Ontology data to list
GO_list <- lapply(GO_list, function(x) unlist(strsplit(x, split='\\`'))) # split single GO terms string into a character vector, one element per term
GO_list <- lapply(GO_list, function(x) substr(x, 1, 10)) # strips away GO descriptions, leaving just ID number for each element
geneID2GO <- GO_list
# make full list of transcript names, geneNames
geneNames <- names(GO_list)

# view first 50 elements
# head(GO_list, 50) 

#splitGOs <- lapply(GO_list, function(x) unlist(strsplit(x, split='\\`'))) # split GO terms string into a character vector, one element per term
#splitGOs_IDonly <- lapply(splitGOs, function(x) substr(x, 1, 10)) # gets just GO ID number for each element


####
# Iterate over all of the topTags result files
#####

# creates the GSEA output file
print("writing output to:")
print(paste(out_dir, out_file, sep=""))
con <- file(paste(out_dir, out_file, sep=""))
sink(con, type = c("output", "message"), split = TRUE)

# iterates over the topTags files, performing GSEA test
for (filename in tTag_output) {
  print("getting file: ")
  print(paste(getwd(), "/", filename, sep=""))
  # Make a list of interesting genes
  # documentation steps 4.4
  # get toptags file from DGE analysis
  DEGs_list <- read.csv(paste(getwd(), "/", filename, sep=""))
  # make a list of interesting transcript IDs based on DEGs_list
  myInterestingGenes2 <- as.character(DEGs_list$transcript_id)
  #head(myInterestingGenes2)
  #str(myInterestingGenes2)
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes2))
  names(geneList) <- geneNames
  # iterate over each of the ontologies
  for (ontol in c("BP", "MF", "CC")) {
    # build the topGOdata object
    GOdata <- new("topGOdata", ontology = ontol, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    # add a decription:
    description(GOdata) <- expDescription
    #description(GOdata)
    #GOdata  
    # generate basic stats of GOdata object:
    GOstats <- termStat(GOdata)
    # to run a test, we need two objects, the topGOdata object already created above, and a groupStats object, which depends on which statistical test you want to run.  See section 6 of the documentation.
    test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
    resultFisher <- getSigGroups(GOdata, test.stat)
    print(resultFisher)  
    allRes <- GenTable(GOdata, classic = resultFisher, topNodes = numNodes)
    print(allRes)
    # Bonferroni correction is simply p-value divided by number of comparisons, in this case, the number of GO terms tested, x in the output where "the algorithm is scoring x nontrivial nodes"
    # build a graph
    #showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = "all")
  }
}

# restores output sink and closes connection to the output file
sink() 
sink(type = c("output", "message"))
rm(con)
