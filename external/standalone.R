# CardioPipeline 2020
# Author: Rohit Suratekar
#
# This is part of the CardioPipeline of the BooleanTRN project
# Standalone R Scrip to perform Differential Expression Analysis.
# Fill in the "sample.yml" file according to your need and run the file
# 
# This stand alone script analyses simple samples

suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("data.table"))

sample <- c("SRR7882015", "SRR7882016", "SRR7882021", "SRR7882022")
condition <- c("24hpf", "24hpf", "72hpf", "72hpf")
reference <- "24hpf"
design <- "~ condition"
count_file <- "/mnt/windows/Enigma/Zebrafish/data/deseq2/mutant/counts/star.counts"
# One can generate count file with the help of bam_to_count.R

colData <- data.frame(sample, condition) # Create the dataframee
df <- read.csv(count_file)
counts <- data.frame(df)
rownames(counts) <- counts$gene_id # Make gene IDs as row names
counts <- counts[sample] # Keep only column equal to our samples in that order

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                      colData = colData,
                                      design = formula(design))

# Relevel the condition by providing information about reference condition
dds$condition <- relevel(dds$condition, ref = reference)
# Start the DE analysis
dds <- DESeq(dds)

# Get results from each genotype combinations
for (c in unique(condition)) {
  if (c != reference) {
    contrast <- c("condition", c, reference)
    res <- results(dds, contrast = contrast)
    name <- paste(c, "vs", reference, sep = "_")
    df <- as.data.frame(res)
    data.table::setDT(df, keep.rownames = "gene_id")
    write.csv(df, file = paste0(name, ".csv"), row.names = FALSE)
  }
}
