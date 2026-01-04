read_merra2_elev.o read_merra2_elev.d : read_merra2_elev.F90
read_merra2_elev.o : LIS_logMod.o
read_merra2_elev.o : LIS_fileIOMod.o
read_merra2_elev.o : LIS_metforcingMod.o
read_merra2_elev.o : LIS_coreMod.o
