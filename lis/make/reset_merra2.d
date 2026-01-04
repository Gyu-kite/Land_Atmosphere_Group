reset_merra2.o reset_merra2.d : reset_merra2.F90
reset_merra2.o : merra2_forcingMod.o
reset_merra2.o : LIS_timeMgrMod.o
reset_merra2.o : LIS_coreMod.o
