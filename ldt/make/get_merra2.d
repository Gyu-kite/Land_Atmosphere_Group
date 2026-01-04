get_merra2.o get_merra2.d : get_merra2.F90
get_merra2.o : merra2_forcingMod.o
get_merra2.o : LDT_timeMgrMod.o
get_merra2.o : LDT_constantsMod.o
get_merra2.o : LDT_logMod.o
get_merra2.o : LDT_metforcingMod.o
get_merra2.o : LDT_coreMod.o
