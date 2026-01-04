read_Precip_climo.o read_Precip_climo.d : read_Precip_climo.F90
read_Precip_climo.o : LDT_misc.h
read_Precip_climo.o : LDT_coreMod.o
read_Precip_climo.o : LDT_historyMod.o
read_Precip_climo.o : LDT_logMod.o
read_Precip_climo.o : LDT_NetCDF_inc.h
