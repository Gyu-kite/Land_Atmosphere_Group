read_Loobos.o read_Loobos.d : read_Loobos.F90
read_Loobos.o : LIS_constantsMod.o
read_Loobos.o : LIS_timeMgrMod.o
read_Loobos.o : LIS_coreMod.o
read_Loobos.o : LIS_metforcingMod.o
read_Loobos.o : Loobos_forcingMod.o
read_Loobos.o : LIS_logMod.o
