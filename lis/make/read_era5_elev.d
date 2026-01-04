read_era5_elev.o read_era5_elev.d : read_era5_elev.F90
read_era5_elev.o : era5_forcingMod.o
read_era5_elev.o : LIS_misc.h
read_era5_elev.o : LIS_coreMod.o
read_era5_elev.o : LIS_logMod.o
read_era5_elev.o : LIS_fileIOMod.o
read_era5_elev.o : LIS_metforcingMod.o
