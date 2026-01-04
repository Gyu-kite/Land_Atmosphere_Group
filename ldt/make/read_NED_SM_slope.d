read_NED_SM_slope.o read_NED_SM_slope.d : read_NED_SM_slope.F90
read_NED_SM_slope.o : LDT_coreMod.o
read_NED_SM_slope.o : LDT_paramTileInputMod.o
read_NED_SM_slope.o : LDT_logMod.o
read_NED_SM_slope.o : LDT_fileIOMod.o
read_NED_SM_slope.o : LDT_gridmappingMod.o
