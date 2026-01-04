read_WRF_AKdom.o read_WRF_AKdom.d : read_WRF_AKdom.F90
read_WRF_AKdom.o : LIS_logMod.o
read_WRF_AKdom.o : LIS_misc.h
read_WRF_AKdom.o : LIS_mpiMod.o
read_WRF_AKdom.o : LIS_coreMod.o
read_WRF_AKdom.o : LIS_constantsMod.o
read_WRF_AKdom.o : LIS_forecastMod.o
read_WRF_AKdom.o : WRF_AKdom_forcingMod.o
read_WRF_AKdom.o : LIS_metforcingMod.o
read_WRF_AKdom.o : LIS_timeMgrMod.o
