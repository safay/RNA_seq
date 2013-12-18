# this script gets data from the DGE .R script and Trinotate output ready for gene set enrichment analysis using the topGO R package
# http://www.bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
# basically this replaces the readMappings() function to take in Trinotate data
# see section 4.3 of the topGO documentation, "custom annotations"

# to install packages:
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("topGO")
# biocLite("Rgraphviz")

# set a character vector as a description of the experiment
expDescription <- "Calineuria thermal stress experiment"
# working directory with your topTags .CSVs from the DGE script
setwd("/Users/Shared/BigCB_Insect_Project/Cali_project/Cali_DGE/")
out_dir <- "/Users/Shared/BigCB_Insect_Project/Cali_project/Cali_DGE/"

#####
# Make full GO list from Trinotate-annotated transcriptome
#####

# get full GO list from transcriptome, a tab-delimited file of Trinotate output, saved in Excel, then processed with the perl one-liner above
annot <- read.delim("/Users/scottfay/annotation/Calineuria/trinotate_annotation_report.txt")
# create a list of GO terms
GO_list <- as.list(as.character(annot[['transcript_id']])) # generates list; element names are transcript IDs
GO_list <- as.list(setNames(as.character(annot$gene_ontology), as.character(annot$transcript_id))) # adds Gene Ontology data to list
GO_list <- lapply(GO_list, function(x) unlist(strsplit(x, split='\\`'))) # split single GO terms string into a character vector, one element per term
GO_list <- lapply(GO_list, function(x) substr(x, 1, 10)) # strips away GO descriptions, leaving just ID number for each element
geneID2GO <- GO_list
# make full list of transcript names, geneNames
geneNames <- names(GO_list)

# view first 50 elements
head(GO_list, 50) 

#splitGOs <- lapply(GO_list, function(x) unlist(strsplit(x, split='\\`'))) # split GO terms string into a character vector, one element per term
#splitGOs_IDonly <- lapply(splitGOs, function(x) substr(x, 1, 10)) # gets just GO ID number for each element


####
# Iterate over all of the topTags result files
#####

# creates a vector of the DGE output files, "tTag_output"
tTag_output <- list.files()
tTag_output <- tTag_output[grepl("topTags.*csv", tTag_output)] # handy little grepl() function!

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
    allRes <- GenTable(GOdata, classic = resultFisher, topNodes = 40)
    print(allRes)
    # Bonferroni correction is simply p-value divided by number of comparisons, in this case, the number of GO terms tested, x in the output where "the algorithm is scoring x nontrivial nodes"
    # build a graph
    #showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = "all")
  }
}



