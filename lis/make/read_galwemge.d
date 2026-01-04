read_galwemge.o read_galwemge.d : read_galwemge.F90
read_galwemge.o : LIS_misc.h
read_galwemge.o : LIS_coreMod.o
read_galwemge.o : LIS_constantsMod.o
read_galwemge.o : galwemge_forcingMod.o
read_galwemge.o : LIS_logMod.o
