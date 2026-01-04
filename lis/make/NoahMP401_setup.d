NoahMP401_setup.o NoahMP401_setup.d : NoahMP401_setup.F90
NoahMP401_setup.o : NoahMP401_lsmMod.o
NoahMP401_setup.o : LIS_coreMod.o
NoahMP401_setup.o : LIS_fileIOMod.o
NoahMP401_setup.o : LIS_misc.h
NoahMP401_setup.o : module_sf_noahmplsm_401.o
NoahMP401_setup.o : LIS_logMod.o
