# CardioPipeline 2020
# Author: Rohit Suratekar
#
# This is part of the CardioPipeline of the BooleanTRN project
# Standalone R Scrip to perform Differential Expression Analysis.
# Fill in the "sample.yml" file according to your need and run the file
# 
# This pipeline currently supports all formats supported by "tximport" package. 
# Hence analysis from BAM files is not supported in this file

library("yaml")
library("tximport")
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("AnnotationDbi"))
suppressPackageStartupMessages(library("data.table"))

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


get_tx2gene <- function(gtf) {
  tx <- GenomicFeatures::makeTxDbFromGFF(file = gtf,
                                         format = "gtf",
                                         organism = "Danio rerio")

  k <- AnnotationDbi::keys(tx, keytype = "GENEID")
  tx_df <- AnnotationDbi::select(tx,
                                 keys = k,
                                 keytype = "GENEID",
                                 columns = "TXNAME")
  return(tx_df[, 2:1])
}

# Get the Transcript to Gene mapping
tx2gene = get_tx2gene(gtf)

# Make the tximport object which counts the all values
txi <- tximport::tximport(files,
                          type = data$method,
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)

# Generate the DDS object
dds <- DESeq2::DESeqDataSetFromTximport(txi = txi,
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
