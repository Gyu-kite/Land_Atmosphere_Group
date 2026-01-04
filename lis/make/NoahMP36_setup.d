NoahMP36_setup.o NoahMP36_setup.d : NoahMP36_setup.F90
NoahMP36_setup.o : module_sf_noahmplsm_36.o
NoahMP36_setup.o : NoahMP36_lsmMod.o
NoahMP36_setup.o : LIS_coreMod.o
NoahMP36_setup.o : LIS_misc.h
NoahMP36_setup.o : LIS_fileIOMod.o
NoahMP36_setup.o : LIS_logMod.o
