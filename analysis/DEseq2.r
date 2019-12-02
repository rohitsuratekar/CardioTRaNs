library(tximport)
library(GenomicFeatures)




gtf_file <- "data/Danio_rerio.GRCz11.98.chr.gtf"
filename <- "data/rnaseq/SRX5113182/salmon/quant.sf"

txdb2 <- makeTxDbFromGFF(file=gtf_file, format="gtf", organism="Danio rerio")

print(tbxdb2)