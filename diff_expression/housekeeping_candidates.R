# finds candidate housekeeping genes 
# based on low coefficient of variation among FPKM values across all libraries
#
# generates a dataframe of mean FPKM, FPKM coefficient of variation

setwd("/Users/Shared/RKC_project/")
load("RKC2_count_and_fpkm_data")

# get annotations file, "trinotate_annotation_report.txt", a Trinotate output excel file exported tab-delimited
transcripts <- read.delim("/Users/scottfay/Dropbox/Red_King_Crab_Project/RKC_transcriptome_and_annotation/trinotate_annotation_report.txt")

# get mean from each transcript
gene.variance.summary <- data.frame(means = apply(all_fpkm, 1, mean))

# get SD from each transcript
gene.variance.summary$SD <- apply(all_fpkm, 1, sd)

# get coefficient of variation from each treatment: ratio of the standard deviation to the mean
gene.variance.summary$CV <- gene.variance.summary$SD / gene.variance.summary$means

# make a new column using the rownames
gene.variance.summary$transcript_id <- rownames(gene.variance.summary)

# merge annotation with summary expression data
merged.transcripts <- merge(transcripts, gene.variance.summary)

# reorder based on coefficient of variation
HK.candidates <- merged.transcripts[order(merged.transcripts$CV),]

# Write a .csv file of the candidate housekeeping genes 
write.csv(head(HK.candidates[,c(1,5,12,13,14)], 100), file="HK_candidates.csv")
