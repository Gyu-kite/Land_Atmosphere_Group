read_nldas2a.o read_nldas2a.d : read_nldas2a.F90
read_nldas2a.o : LDT_misc.h
read_nldas2a.o : LDT_metforcingMod.o
read_nldas2a.o : LDT_spatialDownscalingMod.o
read_nldas2a.o : nldas2_forcingMod.o
read_nldas2a.o : LDT_coreMod.o
read_nldas2a.o : LDT_logMod.o
