reset_gdas.o reset_gdas.d : reset_gdas.F90
reset_gdas.o : LIS_misc.h
reset_gdas.o : gdas_forcingMod.o
reset_gdas.o : LIS_timeMgrMod.o
reset_gdas.o : LIS_coreMod.o
