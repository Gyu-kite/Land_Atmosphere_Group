read_NCEP_GFS_tbot.o read_NCEP_GFS_tbot.d : read_NCEP_GFS_tbot.F90
read_NCEP_GFS_tbot.o : LDT_coreMod.o
read_NCEP_GFS_tbot.o : LDT_logMod.o
read_NCEP_GFS_tbot.o : LDT_gridmappingMod.o
read_NCEP_GFS_tbot.o : Noah_parmsMod.o
