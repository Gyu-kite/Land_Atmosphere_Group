read_gswp2.o read_gswp2.d : read_gswp2.F90
read_gswp2.o : LIS_coreMod.o
read_gswp2.o : LIS_gswpMod.o
read_gswp2.o : gswp2_forcingMod.o
read_gswp2.o : LIS_constantsMod.o
read_gswp2.o : LIS_logMod.o
read_gswp2.o : LIS_metforcingMod.o
read_gswp2.o : LIS_misc.h
