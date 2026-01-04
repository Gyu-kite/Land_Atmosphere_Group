readconfig_WRFoutv2.o readconfig_WRFoutv2.d : readconfig_WRFoutv2.F90
readconfig_WRFoutv2.o : LIS_logMod.o
readconfig_WRFoutv2.o : LIS_coreMod.o
readconfig_WRFoutv2.o : WRFoutv2_forcingMod.o
