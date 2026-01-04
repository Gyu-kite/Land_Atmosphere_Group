readconfig_chirps2.o readconfig_chirps2.d : readconfig_chirps2.F90
readconfig_chirps2.o : LIS_coreMod.o
readconfig_chirps2.o : LIS_logMod.o
readconfig_chirps2.o : chirps2_forcingMod.o
