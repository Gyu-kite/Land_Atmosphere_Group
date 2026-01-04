get_cmorph.o get_cmorph.d : get_cmorph.F90
get_cmorph.o : LDT_timeMgrMod.o
get_cmorph.o : LDT_coreMod.o
get_cmorph.o : LDT_constantsMod.o
get_cmorph.o : cmorph_forcingMod.o
get_cmorph.o : LDT_logMod.o
