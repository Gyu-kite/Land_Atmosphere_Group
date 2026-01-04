read_nldas2b.o read_nldas2b.d : read_nldas2b.F90
read_nldas2b.o : LDT_logMod.o
read_nldas2b.o : nldas2_forcingMod.o
read_nldas2b.o : LDT_coreMod.o
read_nldas2b.o : LDT_metforcingMod.o
read_nldas2b.o : LDT_misc.h
