readcrd_nldas2.o readcrd_nldas2.d : readcrd_nldas2.F90
readcrd_nldas2.o : nldas2_forcingMod.o
readcrd_nldas2.o : LDT_coreMod.o
readcrd_nldas2.o : LDT_logMod.o
