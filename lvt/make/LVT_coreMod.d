LVT_coreMod.o LVT_coreMod.d : LVT_coreMod.F90
LVT_coreMod.o : LVT_PRIV_tileMod.o
LVT_coreMod.o : LVT_mpiMod.o
LVT_coreMod.o : LVT_misc.h
LVT_coreMod.o : LVT_timeMgrMod.o
LVT_coreMod.o : map_utils.o
LVT_coreMod.o : LVT_PRIV_rcMod.o
LVT_coreMod.o : LVT_PRIV_gridMod.o
LVT_coreMod.o : LVT_logMod.o
