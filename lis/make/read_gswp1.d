read_gswp1.o read_gswp1.d : read_gswp1.F90
read_gswp1.o : LIS_logMod.o
read_gswp1.o : LIS_metforcingMod.o
read_gswp1.o : gswp1_forcingMod.o
read_gswp1.o : LIS_coreMod.o
