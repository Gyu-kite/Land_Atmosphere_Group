read_CLM45_lc.o read_CLM45_lc.d : read_CLM45_lc.F90
read_CLM45_lc.o : LDT_coreMod.o
read_CLM45_lc.o : LDT_gridmappingMod.o
read_CLM45_lc.o : LDT_misc.h
read_CLM45_lc.o : CLM45_parmsMod.o
read_CLM45_lc.o : LDT_logMod.o
read_CLM45_lc.o : map_utils.o
