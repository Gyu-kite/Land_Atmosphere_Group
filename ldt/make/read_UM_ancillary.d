read_UM_ancillary.o read_UM_ancillary.d : read_UM_ancillary.F90
read_UM_ancillary.o : LDT_logMod.o
read_UM_ancillary.o : LDT_misc.h
read_UM_ancillary.o : LDT_fileIOMod.o
read_UM_ancillary.o : LDT_paramTileInputMod.o
read_UM_ancillary.o : LDT_gridmappingMod.o
read_UM_ancillary.o : LDT_coreMod.o
