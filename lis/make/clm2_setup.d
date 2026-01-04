clm2_setup.o clm2_setup.d : clm2_setup.F90
clm2_setup.o : clm2_varcon.o
clm2_setup.o : LIS_logMod.o
clm2_setup.o : LIS_timeMgrMod.o
clm2_setup.o : clm2_lsmMod.o
clm2_setup.o : LIS_coreMod.o
