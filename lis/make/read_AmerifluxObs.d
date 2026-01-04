read_AmerifluxObs.o read_AmerifluxObs.d : read_AmerifluxObs.F90
read_AmerifluxObs.o : map_utils.o
read_AmerifluxObs.o : LIS_coreMod.o
read_AmerifluxObs.o : LIS_logMod.o
read_AmerifluxObs.o : LIS_constantsMod.o
read_AmerifluxObs.o : LIS_fileIOMod.o
read_AmerifluxObs.o : LIS_timeMgrMod.o
read_AmerifluxObs.o : AmerifluxobsMod.o
