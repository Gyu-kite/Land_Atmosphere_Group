get_galwem.o get_galwem.d : get_galwem.F90
get_galwem.o : LIS_metforcingMod.o
get_galwem.o : LIS_coreMod.o
get_galwem.o : LIS_logMod.o
get_galwem.o : galwem_forcingMod.o
get_galwem.o : LIS_mpiMod.o
get_galwem.o : LIS_timeMgrMod.o
get_galwem.o : LIS_constantsMod.o
