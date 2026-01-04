readcrd_cmorph.o readcrd_cmorph.d : readcrd_cmorph.F90
readcrd_cmorph.o : LDT_coreMod.o
readcrd_cmorph.o : cmorph_forcingMod.o
readcrd_cmorph.o : LDT_logMod.o
