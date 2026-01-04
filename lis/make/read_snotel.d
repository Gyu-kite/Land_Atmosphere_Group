read_snotel.o read_snotel.d : read_snotel.F90
read_snotel.o : LIS_logMod.o
read_snotel.o : LIS_constantsMod.o
read_snotel.o : LIS_coreMod.o
read_snotel.o : LIS_String_Utility.o
read_snotel.o : LIS_metforcingMod.o
read_snotel.o : map_utils.o
read_snotel.o : snotel_forcingMod.o
