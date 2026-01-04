read_merra2_elev.o read_merra2_elev.d : read_merra2_elev.F90
read_merra2_elev.o : LDT_coreMod.o
read_merra2_elev.o : LDT_metforcingMod.o
read_merra2_elev.o : LDT_logMod.o
read_merra2_elev.o : LDT_misc.h
read_merra2_elev.o : merra2_forcingMod.o
read_merra2_elev.o : LDT_fileIOMod.o
