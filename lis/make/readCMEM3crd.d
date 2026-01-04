readCMEM3crd.o readCMEM3crd.d : readCMEM3crd.F90
readCMEM3crd.o : LIS_misc.h
readCMEM3crd.o : LIS_coreMod.o
readCMEM3crd.o : LIS_logMod.o
readCMEM3crd.o : CMEM3_Mod.o
