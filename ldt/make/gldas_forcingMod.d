gldas_forcingMod.o gldas_forcingMod.d : gldas_forcingMod.F90
gldas_forcingMod.o : LDT_coreMod.o
gldas_forcingMod.o : LDT_logMod.o
gldas_forcingMod.o : LDT_constantsMod.o
gldas_forcingMod.o : LDT_timeMgrMod.o
