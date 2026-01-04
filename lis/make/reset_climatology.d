reset_climatology.o reset_climatology.d : reset_climatology.F90
reset_climatology.o : climatology_forcingMod.o
reset_climatology.o : LIS_misc.h
reset_climatology.o : LIS_coreMod.o
