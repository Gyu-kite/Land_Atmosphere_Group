read_era5_elev.o read_era5_elev.d : read_era5_elev.F90
read_era5_elev.o : LDT_metforcingMod.o
read_era5_elev.o : LDT_logMod.o
read_era5_elev.o : LDT_fileIOMod.o
read_era5_elev.o : LDT_misc.h
read_era5_elev.o : LDT_coreMod.o
read_era5_elev.o : era5_forcingMod.o
