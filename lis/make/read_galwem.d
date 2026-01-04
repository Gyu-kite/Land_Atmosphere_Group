read_galwem.o read_galwem.d : read_galwem.F90
read_galwem.o : LIS_coreMod.o
read_galwem.o : LIS_logMod.o
read_galwem.o : LIS_constantsMod.o
read_galwem.o : galwem_forcingMod.o
read_galwem.o : LIS_misc.h
