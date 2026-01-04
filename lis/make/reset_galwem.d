reset_galwem.o reset_galwem.d : reset_galwem.F90
reset_galwem.o : LIS_misc.h
reset_galwem.o : LIS_coreMod.o
reset_galwem.o : galwem_forcingMod.o
