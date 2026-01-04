read_nldas20_elev.o read_nldas20_elev.d : read_nldas20_elev.F90
read_nldas20_elev.o : nldas20_forcingMod.o
read_nldas20_elev.o : LIS_fileIOMod.o
read_nldas20_elev.o : LIS_metforcingMod.o
read_nldas20_elev.o : LIS_misc.h
read_nldas20_elev.o : LIS_logMod.o
read_nldas20_elev.o : LIS_coreMod.o
