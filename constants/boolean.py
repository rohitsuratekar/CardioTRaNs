#   Copyright (c)  2020, BooleanTrans (previously CardioTrans)
#   Author: Rohit Suratekar
#   Email: rsuratekar [at] iimcb.gov.pl
#   URL: https://github.com/rohitsuratekar/CardioTRaNs
#   Organization: Winata Lab, IIMCB, Warsaw
#
# Constants related to Boolean Modelling


INTERESTED_GENES = ["nkx2.5", "gata4", "hand2", "hey2", "gata2a", "zfpm1",
                    "gata3", "smarca4a", "tbx5a", "mef2ca", "nkx2.3", "foxh1",
                    "gata2b", "sqstm1", "tbx5b", "zfpm2a", "gata5"]

BASE_GENES = ["nkx2.5", "gata4", "tbx5a", "mef2ca", "hand2"]

# Controls from the paper

POSITIVE_CONTROL = ["nkx2.5", "myl7", "myh6"]

NEGATIVE_CONTROL = ["ins", "frem2a", "arr3b", "otx5", "tp63", "col6a1",
                    "neurog1", "zic3", "otx1", "pou4f2", "myog"]

POSITIVE_SET = ['pabpc1a', 'rad21a', 'pdcd10a', 'acta1b', 'desma', 'vdac2',
                'camk2g1', 'unc45b', 'atp1a1a.1', 'cfl1', 'tnni2a.4', 'mat2aa',
                'smarca4a', 'mef2ca', 'rbm24a', 'yap1', 'mef2d', 'gata5',
                'ptch1', 'tpm4a', 'cav1', 'rbfox1l', 'mybpc3', 'cmlc1', 'myl7',
                'fn1a', 'sox9b', 'shha', 'tnnt2a', 'klf2b']

NEGATIVE_SET = ['cyp2n13', 'tagapa', 'nr5a5', 'trpm6', 'c6ast3', 'slc35a3b',
                'il12a', 'poln', 'tdo2b', 'si:ch211-207i1.2', 'si:ch73-7i4.2',
                'arl3a']

HOUSE_KEEPING = ['actb1', 'actb2', 'b2m', 'ppiaa', "rpl32", "h2az2a",
                 "gapdh", "myl7"]

HK_STAR_WINATA = ['si:ch211-248a14.8', 'BX530407.4', 'gsx2', 'zgc:101765',
                  'sptan1', 'dip2a', 'nckap5l', 'cbfb']

HK_STAR_HILLS = ['pbx4', 'zeb2a', 'ctbp1', 'cirbpa', 'marcksa']
HK_SALMON_WINATA = ['ccdc17', 'defbl1', 'prlh', 'amigo3', 'tmem63c']
