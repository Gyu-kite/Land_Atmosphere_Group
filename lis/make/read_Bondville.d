read_Bondville.o read_Bondville.d : read_Bondville.F90
read_Bondville.o : LIS_coreMod.o
read_Bondville.o : LIS_timeMgrMod.o
read_Bondville.o : LIS_constantsMod.o
read_Bondville.o : LIS_logMod.o
read_Bondville.o : LIS_metforcingMod.o
read_Bondville.o : Bondville_forcingMod.o
