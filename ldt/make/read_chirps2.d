read_chirps2.o read_chirps2.d : read_chirps2.F90
read_chirps2.o : LDT_timeMgrMod.o
read_chirps2.o : LDT_logMod.o
read_chirps2.o : chirps2_forcingMod.o
read_chirps2.o : LDT_misc.h
read_chirps2.o : LDT_coreMod.o
read_chirps2.o : LDT_metforcingMod.o
