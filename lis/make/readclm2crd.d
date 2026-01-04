readclm2crd.o readclm2crd.d : readclm2crd.F90
readclm2crd.o : LIS_logMod.o
readclm2crd.o : LIS_timeMgrMod.o
readclm2crd.o : clm2_lsmMod.o
readclm2crd.o : LIS_coreMod.o
