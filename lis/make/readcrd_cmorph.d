readcrd_cmorph.o readcrd_cmorph.d : readcrd_cmorph.F90
readcrd_cmorph.o : LIS_logMod.o
readcrd_cmorph.o : LIS_coreMod.o
readcrd_cmorph.o : cmorph_forcingMod.o
