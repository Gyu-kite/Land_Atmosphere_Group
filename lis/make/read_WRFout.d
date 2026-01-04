read_WRFout.o read_WRFout.d : read_WRFout.F90
read_WRFout.o : WRFout_forcingMod.o
read_WRFout.o : LIS_misc.h
read_WRFout.o : LIS_logMod.o
read_WRFout.o : LIS_metforcingMod.o
read_WRFout.o : LIS_coreMod.o
