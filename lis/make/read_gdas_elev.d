read_gdas_elev.o read_gdas_elev.d : read_gdas_elev.F90
read_gdas_elev.o : LIS_metforcingMod.o
read_gdas_elev.o : LIS_fileIOMod.o
read_gdas_elev.o : LIS_coreMod.o
read_gdas_elev.o : gdas_forcingMod.o
read_gdas_elev.o : LIS_logMod.o
