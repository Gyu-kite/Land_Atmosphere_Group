read_cmap.o read_cmap.d : read_cmap.F90
read_cmap.o : LDT_coreMod.o
read_cmap.o : LDT_logMod.o
read_cmap.o : LDT_metforcingMod.o
read_cmap.o : LDT_misc.h
read_cmap.o : cmap_forcingMod.o
