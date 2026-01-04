read_mogrepsg.o read_mogrepsg.d : read_mogrepsg.F90
read_mogrepsg.o : LIS_coreMod.o
read_mogrepsg.o : LIS_constantsMod.o
read_mogrepsg.o : mogrepsg_forcingMod.o
read_mogrepsg.o : LIS_logMod.o
read_mogrepsg.o : LIS_misc.h
