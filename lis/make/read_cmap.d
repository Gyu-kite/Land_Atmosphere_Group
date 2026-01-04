read_cmap.o read_cmap.d : read_cmap.F90
read_cmap.o : cmap_forcingMod.o
read_cmap.o : LIS_misc.h
read_cmap.o : LIS_coreMod.o
read_cmap.o : LIS_logMod.o
read_cmap.o : LIS_metforcingMod.o
