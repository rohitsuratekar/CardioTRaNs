# CardioPipeline 2020
# Author: Rohit Suratekar
#
# This is part of the CardioPipeline of the BooleanTRN project
# Standalone R Scrip to perform Differential Expression Analysis.
# Fill in the "sample.yml" file according to your need and run the file
# 
# This pipeline currently supports RNA-seq DE analysis on the BAM files (usually generated with STAR)

library("yaml")
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("data.table"))
library("Rsubread")

# Read the sample file and create DataFrame with sample names and conditions
data <- yaml::read_yaml("samples.yml")
gtf <- data$gtf
sample <- c()
genotype <- c()
files <- c()
for (key in names(data$sample)) {
  sample <- c(sample, key)
  gn <- data$sample[[key]][["genotype"]]
  fl <- data$samples[[key]][["input"]]
  genotype <- c(genotype, gn)
  files <- c(files, fl)
}
names(files) <- sample  # Give column names are respective sample names
colData <- data.frame(sample, genotype) # Create the dataframee

# Your count file which will have gene names on rows and column names are sample IDs.
# You can generate this from BAM files with script 'bam_to_count.R'
count_file <- "star.counts"
df <- read.csv(count_file)
counts <- data.frame(df)
rownames(counts) <- counts$gene_id # Make gene IDs as row names
counts$gene_id <- NULL # Remove the column
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                      colData = colData,
                                      design = formula(data$design))
# Relevel the condition by providing information about reference condition
dds$genotype <- relevel(dds$genotype, ref = data$reference)
# Start the DE analysis
dds <- DESeq(dds)

# Get results from each genotype combinations
for (g in unique(genotype)) {
  if (g != data$reference) {
    contrast <- c("genotype", g, data$reference)
    res <- results(dds, contrast = contrast)
    name <- paste(data$method
      , g, "vs", data$reference, sep = "_")
    df <- as.data.frame(res)
    data.table::setDT(df, keep.rownames = "gene_id")
    write.csv(df, file = paste0(name, ".csv"), row.names = FALSE)
  }
}
