LDT_coreMod.o LDT_coreMod.d : LDT_coreMod.F90
LDT_coreMod.o : LDT_misc.h
LDT_coreMod.o : map_utils.o
LDT_coreMod.o : LDT_mpiMod.o
LDT_coreMod.o : LDT_PRIV_gridMod.o
LDT_coreMod.o : LDT_timeMgrMod.o
LDT_coreMod.o : LDT_logMod.o
LDT_coreMod.o : LDT_PRIV_rcMod.o
LDT_coreMod.o : LDT_PRIV_tileMod.o
