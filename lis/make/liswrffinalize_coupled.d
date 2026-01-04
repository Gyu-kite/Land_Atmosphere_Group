liswrffinalize_coupled.o liswrffinalize_coupled.d : liswrffinalize_coupled.F90
liswrffinalize_coupled.o : LIS_timeMgrMod.o
liswrffinalize_coupled.o : LIS_logMod.o
liswrffinalize_coupled.o : LIS_irrigationMod.o
liswrffinalize_coupled.o : LIS_coreMod.o
liswrffinalize_coupled.o : LIS_surfaceModelMod.o
liswrffinalize_coupled.o : LISWRFGridCompMod.o
