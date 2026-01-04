read_cmorph.o read_cmorph.d : read_cmorph.F90
read_cmorph.o : LDT_metforcingMod.o
read_cmorph.o : LDT_coreMod.o
read_cmorph.o : LDT_constantsMod.o
read_cmorph.o : LDT_logMod.o
read_cmorph.o : cmorph_forcingMod.o
