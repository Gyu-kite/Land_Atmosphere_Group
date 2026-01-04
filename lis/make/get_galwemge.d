get_galwemge.o get_galwemge.d : get_galwemge.F90
get_galwemge.o : galwemge_forcingMod.o
get_galwemge.o : LIS_timeMgrMod.o
get_galwemge.o : LIS_coreMod.o
get_galwemge.o : LIS_constantsMod.o
get_galwemge.o : LIS_logMod.o
get_galwemge.o : LIS_metforcingMod.o
