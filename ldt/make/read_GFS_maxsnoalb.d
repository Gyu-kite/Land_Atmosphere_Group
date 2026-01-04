read_GFS_maxsnoalb.o read_GFS_maxsnoalb.d : read_GFS_maxsnoalb.F90
read_GFS_maxsnoalb.o : LDT_fileIOMod.o
read_GFS_maxsnoalb.o : LDT_albedoMod.o
read_GFS_maxsnoalb.o : LDT_logMod.o
read_GFS_maxsnoalb.o : LDT_gridmappingMod.o
read_GFS_maxsnoalb.o : LDT_coreMod.o
