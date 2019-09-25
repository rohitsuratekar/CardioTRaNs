#  Copyright (c) 2019.
#  BooleanTRN (previously CardioTraNs)
#  Author: Rohit Suratekar
#  Website: https://github.com/rohitsuratekar/CardioTRaNs
#
# Constants related to NGS Analysis


# Method outputs
OUTPUT_STRING_TIE = "StringTie"
OUTPUT_STAR = "STAR"
OUTPUT_SALMON = "Salmon"
OUTPUT_SALMON_LOG = "Salmon_log"

# Referenced ID for lab databases
BIO_PROJECT_WINATA_LAB = "PRJNA492280"
BIO_PROJECT_COUGHLIN_LAB = "PRJNA474426"

# RNA seq related
R_SEQ_FILE_SRA_ID = "sra_id"
R_SEQ_FILE_BIO_PROJECT = "bioproject"
R_SEQ_FILE_TIME = "time"
R_SEQ_FILE_GENOTYPE = "genotype"
R_SEQ_FILE_TISSUE = "tissue"
R_SEQ_FILE_INSTRUMENT = "instrument"
R_SEQ_FILE_LAYOUT = "layout"
R_SEQ_FILE_SELECTION = "selection"
R_SEQ_FILE_REFERENCE = "reference"
R_SEQ_FILE_REF_DETAILS = "ref_details"

# STAR output

STAR_JOB_START = "Started job on"
STAR_MAPPING_START = "Started mapping on"
STAR_FINISHED_ON = "Finished on"
STAR_MAPPING_SPEED = "Mapping speed, Million of reads per hour"
STAR_INPUT_READS = "Number of input reads"
STAR_AVG_READ_LENGTH = "Average input read length"
STAR_UNIQUELY_MAPPED_READS = "Uniquely mapped reads number"
STAR_UNIQUELY_MAPPED_READ_PERCENTAGE = "Uniquely mapped reads %"
STAR_AVG_MAPPED_LENGTH = "Average mapped length"
STAR_TOTAL_SPLICES = "Number of splices: Total"
STAR_ANNOTATED_SPLICES = "Number of splices: Annotated (sjdb)"
STAR_SPLICES_GT_AG = "Number of splices: GT/AG"
STAR_SPLICES_GC_AG = "Number of splices: GC/AG"
STAR_SPLICES_AT_AC = "Number of splices: AT/AC"
STAR_SPLICES_NON_CANONICAL = "Number of splices: Non-canonical"
STAR_MISMATCH_RATE_PERCENT = "Mismatch rate per base, %"
STAR_DELETION_RATE = "Deletion rate per base"
STAR_AVG_DELETION_LENGTH = "Deletion average length"
STAR_INSERTION_RATE = "Insertion rate per base"
STAR_AVG_INSERTION_LENGTH = "Insertion average length"
STAR_MULTI_MAPPED_READS = "Number of reads mapped to multiple loci"
STAR_MULTI_MAPPED_PERCENT = "% of reads mapped to multiple loci"
STAR_TOO_MANY_LOCI_MAPPING = "Number of reads mapped to too many loci"
STAR_TOO_MANY_LOCI_MAPPING_PERCENT = "% of reads mapped to too many loci"
STAR_MISMATCH_UNMAPPED_READS = "Number of reads unmapped: too many mismatches"
STAR_MISMATCH_UNMAPPED_READS_PERCENT = "% of reads unmapped: too many mismatches"
STAR_SHORT_UNMAPPED_READS = "Number of reads unmapped: too short"
STAR_SHORT_UNMAPPED_READS_PERCENT = "% of reads unmapped: too short"
STAR_OTHER_UNMAPPED_READS = "Number of reads unmapped: other"
STAR_OTHER_UNMAPPED_READS_PERCENT = "% of reads unmapped: other"
STAR_CHIMERIC_READS = "Number of chimeric reads"
STAR_CHIMERIC_READS_PERCENT = "% of chimeric reads"

# Salmon output file columns

SALMON_NAME = "Name"
SALMON_LENGTH = "Length"
SALMON_EFFECTIVE_LENGTH = "EffectiveLength"
SALMON_TPM = "TPM"
SALMON_NUM_READS = "NumReads"
SALMON_GENE_STABLE_ID = "Gene stable ID"
SALMON_GENE_NAME = "Gene name"
SALMON_ZFIN_ID = "ZFIN ID"
SALMON_TRANSCRIPT_SOURCE = "Source (transcript)"

# Salmon Log

SAL_LOG_FRAG_DIST_LENGTH = "frag_dist_length"
SAL_LOG_VALID_TARGETS = "num_valid_targets"
SAL_LOG_DECOY_TARGETS = "num_decoy_targets"
SAL_LOG_PROCESSED_READS = "num_processed"
SAL_LOG_MAPPED_READS = "num_mapped"
SAL_LOG_MAPPING_PERCENT = "percent_mapped"
SAL_LOG_START_TIME = "start_time"
SAL_LOG_END_TIME = "end_time"

# StringTie constants

STRING_GENE_ID = "Gene ID"
STRING_GENE_NAME = "Gene Name"
STRING_REFERENCE = "Reference"
STRING_STRAND = "Strand"
STRING_START = "Start"
STRING_END = "End"
STRING_COVERAGE = "Coverage"
STRING_FPKM = "FPKM"
STRING_TPM = "TPM"
