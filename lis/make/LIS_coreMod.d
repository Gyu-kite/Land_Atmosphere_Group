LIS_coreMod.o LIS_coreMod.d : LIS_coreMod.F90
LIS_coreMod.o : LIS_PRIV_gridMod.o
LIS_coreMod.o : LIS_timeMgrMod.o
LIS_coreMod.o : LIS_logMod.o
LIS_coreMod.o : LIS_misc.h
LIS_coreMod.o : LIS_PRIV_rcMod.o
LIS_coreMod.o : LIS_PRIV_tileMod.o
LIS_coreMod.o : map_utils.o
LIS_coreMod.o : LIS_mpiMod.o
