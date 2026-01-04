read_imerg.o read_imerg.d : read_imerg.F90
read_imerg.o : LIS_misc.h
read_imerg.o : LIS_metforcingMod.o
read_imerg.o : LIS_coreMod.o
read_imerg.o : LIS_constantsMod.o
read_imerg.o : LIS_logMod.o
read_imerg.o : imerg_forcingMod.o
