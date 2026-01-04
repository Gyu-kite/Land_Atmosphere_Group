reset_gddp.o reset_gddp.d : reset_gddp.F90
reset_gddp.o : LIS_coreMod.o
reset_gddp.o : gddp_forcingMod.o
reset_gddp.o : LIS_misc.h
