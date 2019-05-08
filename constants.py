"""
CardioTRaNs 2019
Author: Rohit Suratekar

All constants of this project (including database constants)
"""

DATABASE = "CardioTrans.hdf5"  # Name of the main HDF5 database
DATABASE_VERSION = 1
DATABASE_AUTHOR = "Rohit Suratekar"

DATA_FOLDER = "data/"  # Database Folder

# All database files
FILE_ANATOMY_ITEMS = "anatomy_item_2019.04.08.txt"
FILE_ANATOMY_RELATIONSHIP = "anatomy_relationship_2019.04.08.txt"
FILE_STAGE_ONTOLOGY = "stage_ontology_2019.04.08.txt"
FILE_WILD_TYPE_EXPRESSION = "wildtype-expression_fish_2019.04.08.txt"
FILE_DE_SEQ_GFP = "Supplement_Table_1.xlsx"

# Column Names for Expression Files
COL_EXP_0_GENE_ID = "gene_id"
COL_EXP_1_GENE_SYMBOL = "gene_symbol"
COL_EXP_2_FISH_NAME = "fish_name"
COL_EXP_3_SUPER_STR_ID = "super_structure_id"
COL_EXP_4_SUPER_STR_NAME = "super_structure_name"
COL_EXP_5_SUB_STR_ID = "sub_structure_id"
COL_EXP_6_SUB_STR_NAME = "sub_structure_name"
COL_EXP_7_START_STAGE = "start_stage"
COL_EXP_8_END_STAGE = "end_stage"
COL_EXP_9_ASSAY = "assay"
COL_EXP_10_ASSAY_ID = "assay_id"
COL_EXP_11_PUBLICATION_ID = "publication_id"
COL_EXP_12_PROBE_ID = "probe_id"
COL_EXP_13_ANTIBODY_ID = "antibody_id"
COL_EXP_14_FISH_ID = "fish_id"

COL_EXP_ALL = [COL_EXP_0_GENE_ID,
               COL_EXP_1_GENE_SYMBOL,
               COL_EXP_2_FISH_NAME,
               COL_EXP_3_SUPER_STR_ID,
               COL_EXP_4_SUPER_STR_NAME,
               COL_EXP_5_SUB_STR_ID,
               COL_EXP_6_SUB_STR_NAME,
               COL_EXP_7_START_STAGE,
               COL_EXP_8_END_STAGE,
               COL_EXP_9_ASSAY,
               COL_EXP_10_ASSAY_ID,
               COL_EXP_11_PUBLICATION_ID,
               COL_EXP_12_PROBE_ID,
               COL_EXP_13_ANTIBODY_ID,
               COL_EXP_14_FISH_ID]

# Columns in Pawlak et. al. 2019 Sup Table 1

COL_GR_S1_GENE_ID = "geneID"
COL_GR_S1_GENE_SYMBOL = "ZFIN_ID"  # In original document, Gene Symbols are
# wrongly labeled as ZFIN ID
COL_GR_S1_24_LOG2FC = "RNAseq_24_plus_minus_log2FC"
COL_GR_S1_48_LOG2FC = "RNAseq_48_plus_minus_log2FC"
COL_GR_S1_72_LOG2FC = "RNAseq_72_plus_minus_log2FC"

# Database Groups
DB_DATA = "zfin_data"
DB_ANATOMY = "anatomy"
DB_ANATOMY_RELATIONSHIPS = "anatomy_relationships"
DB_STAGE_ONTOLOGY = "stage_ontology"
DB_WILD_TYPE_EXP = "wild_type_expression"

# Anatomy constants
DB_ANATOMY_ID = 0
DB_ANATOMY_NAME = 1
DB_ANATOMY_START_STAGE = 2
DB_ANATOMY_END_STAGE = 3
DB_ANATOMY_PARENT_ID = 4

# Stage Constants
DB_STAGE_ID = 0
DB_STAGE_OBO_ID = 1
DB_STAGE_NAME = 2
DB_STAGE_START = 3
DB_STAGE_END = 4

# Expression Constants

DB_EXP_ID = 0
DB_EXP_SYMBOL = 1
DB_EXP_FISH_NAME = 2
DB_EXP_SUPER_STR_ID = 3
DB_EXP_SUB_STR_ID = 4
DB_EXP_START_STAGE_ID = 5
DB_EXP_END_STAGE_ID = 6
DB_EXP_ASSAY = 7
DB_EXP_ASSAY_ID = 8
DB_EXP_PUBLICATION_ID = 9
DB_EXP_PROBE_ID = 10
DB_EXP_ANTIBODY_ID = 11
DB_EXP_FISH_ID = 12

# Meta-data attributes
DB_ATTRS_VERSION = "version"
DB_ATTRS_CREATION = "creation"
DB_ATTRS_AUTHOR = "author"
