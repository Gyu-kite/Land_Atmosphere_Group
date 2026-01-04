read_era5.o read_era5.d : read_era5.F90
read_era5.o : LDT_coreMod.o
read_era5.o : LDT_FORC_AttributesMod.o
read_era5.o : LDT_misc.h
read_era5.o : LDT_metforcingMod.o
read_era5.o : era5_forcingMod.o
read_era5.o : LDT_logMod.o
