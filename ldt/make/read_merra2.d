read_merra2.o read_merra2.d : read_merra2.F90
read_merra2.o : LDT_coreMod.o
read_merra2.o : LDT_FORC_AttributesMod.o
read_merra2.o : LDT_logMod.o
read_merra2.o : LDT_misc.h
read_merra2.o : LDT_metforcingMod.o
read_merra2.o : merra2_forcingMod.o
