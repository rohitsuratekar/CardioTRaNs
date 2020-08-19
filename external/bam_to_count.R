# CardioPipeline 2020
# Author: Rohit Suratekar
#
# This is part of the CardioPipeline of the BooleanTRN project
# Standalone R Scrip to perform Differential Expression Analysis.
# Fill in the "sample.yml" file according to your need and run the file
# 
# Simple script generate counts from the BAM files

library("yaml")
suppressPackageStartupMessages(library("data.table"))
library("Rsubread")

# Read the sample file and create DataFrame with sample names and conditions
data <- yaml::read_yaml("samples.yml")
gtf <- data$gtf
sample <- c()
genotype <- c()
files <- c()
for(key in names(data$sample)){
  sample <- c(sample, key)
  gn <- data$sample[[key]][["genotype"]]
  fl <- data$samples[[key]][["input"]]
  genotype <- c(genotype, gn)
  files <- c(files, fl)
}
names(files) <- sample  # Give column names are respective sample names
colData <- data.frame(sample,genotype) # Create the dataframee

# Change following accordingly
paried_end <- TRUE
threads <- 8

fc <- Rsubread::featureCounts(
  files = files,
  annot.ext = data$gtf,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = paried_end,
  nthreads = threads)

colnames(fc$counts) <- sample
df <- data.frame(fc$counts)
data.table::setDT(df, keep.rownames = "gene_id")
write.csv(df, file= "star.counts",row.names = FALSE)
write.csv(fc$stat, file = "star.counts.stat", row.names = FALSE)
