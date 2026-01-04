read_agrradps.o read_agrradps.d : read_agrradps.F90
read_agrradps.o : LIS_constantsMod.o
read_agrradps.o : LIS_logMod.o
read_agrradps.o : LIS_coreMod.o
read_agrradps.o : LIS_FORC_AttributesMod.o
read_agrradps.o : LIS_metforcingMod.o
read_agrradps.o : agrradps_forcingMod.o
