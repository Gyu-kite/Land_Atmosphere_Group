read_SSMISNWDsnow.o read_SSMISNWDsnow.d : read_SSMISNWDsnow.F90
read_SSMISNWDsnow.o : LIS_timeMgrMod.o
read_SSMISNWDsnow.o : LIS_coreMod.o
read_SSMISNWDsnow.o : LIS_constantsMod.o
read_SSMISNWDsnow.o : LIS_mpiMod.o
read_SSMISNWDsnow.o : LIS_DAobservationsMod.o
read_SSMISNWDsnow.o : LIS_misc.h
read_SSMISNWDsnow.o : LIS_pluginIndices.o
read_SSMISNWDsnow.o : SSMISNWDsnow_Mod.o
read_SSMISNWDsnow.o : LIS_logMod.o
