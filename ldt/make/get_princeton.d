get_princeton.o get_princeton.d : get_princeton.F90
get_princeton.o : LDT_timeMgrMod.o
get_princeton.o : LDT_logMod.o
get_princeton.o : LDT_metforcingMod.o
get_princeton.o : princeton_forcingMod.o
get_princeton.o : LDT_coreMod.o
