read_NCEP_tbot.o read_NCEP_tbot.d : read_NCEP_tbot.F90
read_NCEP_tbot.o : LDT_fileIOMod.o
read_NCEP_tbot.o : LDT_logMod.o
read_NCEP_tbot.o : LDT_coreMod.o
read_NCEP_tbot.o : Noah_parmsMod.o
