get_nldas2.o get_nldas2.d : get_nldas2.F90
get_nldas2.o : nldas2_forcingMod.o
get_nldas2.o : LDT_logMod.o
get_nldas2.o : LDT_constantsMod.o
get_nldas2.o : LDT_metforcingMod.o
get_nldas2.o : LDT_coreMod.o
get_nldas2.o : LDT_timeMgrMod.o
