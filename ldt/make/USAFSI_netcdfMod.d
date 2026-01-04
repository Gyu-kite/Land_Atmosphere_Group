USAFSI_netcdfMod.o USAFSI_netcdfMod.d : USAFSI_netcdfMod.F90
USAFSI_netcdfMod.o : USAFSI_paramsMod.o
USAFSI_netcdfMod.o : LDT_coreMod.o
USAFSI_netcdfMod.o : LDT_logMod.o
USAFSI_netcdfMod.o : LDT_usafsiMod.o
USAFSI_netcdfMod.o : USAFSI_arraysMod.o
USAFSI_netcdfMod.o : LDT_NetCDF_inc.h
USAFSI_netcdfMod.o : LDT_constantsMod.o
USAFSI_netcdfMod.o : LDT_misc.h
