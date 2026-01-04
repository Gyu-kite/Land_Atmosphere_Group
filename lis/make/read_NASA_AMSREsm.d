read_NASA_AMSREsm.o read_NASA_AMSREsm.d : read_NASA_AMSREsm.F90
read_NASA_AMSREsm.o : LIS_pluginIndices.o
read_NASA_AMSREsm.o : LIS_coreMod.o
read_NASA_AMSREsm.o : LIS_timeMgrMod.o
read_NASA_AMSREsm.o : LIS_DAobservationsMod.o
read_NASA_AMSREsm.o : LIS_misc.h
read_NASA_AMSREsm.o : NASA_AMSREsm_Mod.o
read_NASA_AMSREsm.o : LIS_logMod.o
read_NASA_AMSREsm.o : LIS_mpiMod.o
read_NASA_AMSREsm.o : LIS_constantsMod.o
