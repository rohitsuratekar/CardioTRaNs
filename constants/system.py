#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 19/07/19, 2:12 PM
#
#  Copyright (c) 2019.
#
#  All system constants will go in this file


DATA_FOLDER = "data/"
PLOT_FOLDER = "plots/"

GENOTYPE_WT = "wt"
GENOTYPE_GATA5 = "gata5"

FILE_MANUAL = "manual.csv"
FILE_EXPRESSION_ATLAS = "E-ERAD-475-query-results.tpms.tsv"

FILE_ZFIN_EXPRESSION = "zfin/wildtype-expression_fish_2019.07.27.txt"
FILE_ZFIN_ANATOMY_ITEMS = "zfin/anatomy_item_2019.07.27.txt"
FILE_ZFIN_ANATOMY_RELATIONSHIP = "zfin/anatomy_relationship_2019.07.27.txt"
FILE_ZFIN_STAGE_ONTOLOGY = "zfin/stage_ontology_2019.07.27.txt"
FILE_ZFIN_XPAT_STAGE = "zfin/xpat_stage_anatomy_2019.07.27.txt"

# Wild type sequences
FILE_RNA_SEQ_20H_1 = "rna_seq/Z_20h_wt_SRX4157236.tsv"

FILE_RNA_SEQ_24H_1 = "rna_seq/Z_24h_wt_SRX4720625.tsv"
FILE_RNA_SEQ_24H_2 = "rna_seq/Z_24h_wt_SRX4720626.tsv"
FILE_RNA_SEQ_48H_1 = "rna_seq/Z_48h_wt_SRX4720628.tsv"
FILE_RNA_SEQ_48H_2 = "rna_seq/Z_48h_wt_SRX4720629.tsv"
FILE_RNA_SEQ_72H_1 = "rna_seq/Z_72h_wt_SRX4720631.tsv"
FILE_RNA_SEQ_72H_2 = "rna_seq/Z_72h_wt_SRX4720632.tsv"

# Mutant sequences
FILE_MT_GATA5_72H_1 = "rna_seq/Z_72h_gata5_SRX4720634.tsv"
FILE_MT_GATA5_72H_2 = "rna_seq/Z_72h_gata5_SRX4720635.tsv"

# Salmon Files
FILE_SALMON_24H_1 = "salmon/S_24h_WT_SRX4720626.tsv"
FILE_SALMON_48H_1 = "salmon/S_48h_WT_SRX4720628.tsv"
FILE_SALMON_48H_2 = "salmon/S_48h_WT_SRX4720629.tsv"
FILE_SALMON_72H_1 = "salmon/S_72h_WT_SRX4720631.tsv"
FILE_SALMON_72H_2 = "salmon/S_72h_WT_SRX4720632.tsv"

FILE_SALMON_GATA5_72H_1 = "salmon/S_72h_GATA5_SRX4720634.tsv"
FILE_SALMON_GATA5_72H_2 = "salmon/S_72h_GATA5_SRX4720635.tsv"
FILE_SALMON_TBX5_72H_1 = "salmon/S_72h_TBX5_SRX4720637.tsv"

# Mapping log folder from STAR alignment (e.g. SRX4157236_Log.final.out)
# This folder will contain all log files from the STAR alignment
MAPPING_FOLDER = "data/mapping/"
MAPPING_FILE_SUFFIX = "_Log.final.out"

FILE_BOOL_RNA_SEQ = "rna_seq/average_bool_"
