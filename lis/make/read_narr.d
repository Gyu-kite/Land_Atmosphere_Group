read_narr.o read_narr.d : read_narr.F90
read_narr.o : narr_forcingMod.o
read_narr.o : LIS_coreMod.o
read_narr.o : LIS_logMod.o
