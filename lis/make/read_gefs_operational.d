read_gefs_operational.o read_gefs_operational.d : read_gefs_operational.F90
read_gefs_operational.o : gefs_forcingMod.o
read_gefs_operational.o : LIS_coreMod.o
read_gefs_operational.o : LIS_logMod.o
read_gefs_operational.o : LIS_misc.h
read_gefs_operational.o : LIS_metforcingMod.o
