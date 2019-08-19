#  Project: CardioTrans
#  Author: Rohit Suratekar
#  Created On: 29/07/19, 2:10 PM
#
#  Copyright (c) 2019.
#
# All Boolean analysis related constants will go here


# In RNA-seq data, below this cut-off genes will be considered not expressed
TPM_CUT_OFF = 10
FPKM_CUT_OFF = 10

# All genes interested in the boolean model. Based on this list, new small
# file will be generated which can be easily read multiple times

INTERESTED_GENES = ["nkx2.5", "gata4", "hand2", "hey2", "gata2a", "zfpm1",
                    "gata3", "smarca4a", "tbx5a", "mef2ca", "nkx2.3", "foxh1",
                    "gata2b", "sqstm1", "tbx5b", "zfpm2a", "gata5"]

BASE_GENES = ["nkx2.5", "gata4", "tbx5a", "mef2ca", "hand2"]

BOOL_OBJ_COL = "BoolObj"

CONTROL_GENES = ["neurog1", "ins", "arr3b"]
