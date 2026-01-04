read_GTOPO30_GFS_elev.o read_GTOPO30_GFS_elev.d : read_GTOPO30_GFS_elev.F90
read_GTOPO30_GFS_elev.o : LDT_logMod.o
read_GTOPO30_GFS_elev.o : LDT_misc.h
read_GTOPO30_GFS_elev.o : LDT_fileIOMod.o
read_GTOPO30_GFS_elev.o : LDT_gridmappingMod.o
read_GTOPO30_GFS_elev.o : LDT_coreMod.o
