#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#

STRING_COLUMNS = ['protein1', 'protein2', 'neighborhood',
                  'neighborhood_transferred', 'fusion', 'cooccurence',
                  'homology', 'coexpression', 'coexpression_transferred',
                  'experiments', 'experiments_transferred', 'database',
                  'database_transferred', 'textmining',
                  'textmining_transferred', 'combined_score']

# StringTie

STRINGTIE_GENE_ID = "Gene ID"
STRINGTIE_NAME = "Gene Name"
STRINGTIE_REFERENCE = "Reference"
STRINGTIE_STRAND = "Strand"
STRINGTIE_START = "Start"
STRINGTIE_END = "End"
STRINGTIE_COVERAGE = "Coverage"
STRINGTIE_FPKM = "FPKM"
STRINGTIE_TPM = "TPM"

# Salmon

SALMON_NAME = "Name"
SALMON_LENGTH = "Length"
SALMON_EFFECTIVE_LENGTH = "EffectiveLength"
SALMON_TPM = "TPM"
SALMON_NUM_READS = "NumReads"

# Kallisto
KALLISTO_TARGET_ID = "target_id"
KALLISTO_LENGTH = "length"
KALLISTO_EFFECTIVE_LENGTH = "eff_length"
KALLISTO_EST_COUNTS = "est_counts"
KALLISTO_TPM = "tpm"

# DESEq2

DESEQ2_GENE_ID = "gene_id"
DESEQ2_BASE_MEAN = "baseMean"
DESEQ2_LOG2_CHANGE = "log2FoldChange"
DESEQ2_LFCSE = "lfcSE"
DESEQ2_STAT = "stat"
DESEQ2_PVALUE = "pvalue"
DESEQ2_PADJ = "padj"

# Biomart
BIOMART_GENE_ID = 'Gene stable ID'
BIOMART_GENE_ID_VERSION = 'Gene stable ID version'
BIOMART_TRANSCRIPT_ID = 'Transcript stable ID'
BIOMART_TRANSCRIPT_ID_VERSION = 'Transcript stable ID version'
BIOMART_GENE_NAME = 'Gene name'
BIOMART_ZFIN_ID = 'ZFIN ID'
BIOMART_ZFIN_SYMBOL = 'ZFIN symbol'
