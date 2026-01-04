readconfig_chirps2.o readconfig_chirps2.d : readconfig_chirps2.F90
readconfig_chirps2.o : chirps2_forcingMod.o
readconfig_chirps2.o : LDT_logMod.o
readconfig_chirps2.o : LDT_coreMod.o
