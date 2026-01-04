get_climatology.o get_climatology.d : get_climatology.F90
get_climatology.o : climatology_VariablesMod.o
get_climatology.o : climatology_forcingMod.o
get_climatology.o : LIS_coreMod.o
get_climatology.o : LIS_logMod.o
get_climatology.o : LIS_timeMgrMod.o
get_climatology.o : LIS_constantsMod.o
