read_GFS_slopetype.o read_GFS_slopetype.d : read_GFS_slopetype.F90
read_GFS_slopetype.o : LDT_fileIOMod.o
read_GFS_slopetype.o : LDT_logMod.o
read_GFS_slopetype.o : Noah_parmsMod.o
read_GFS_slopetype.o : LDT_gridmappingMod.o
read_GFS_slopetype.o : LDT_coreMod.o
