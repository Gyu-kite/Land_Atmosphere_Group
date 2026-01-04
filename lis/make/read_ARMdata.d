read_ARMdata.o read_ARMdata.d : read_ARMdata.F90
read_ARMdata.o : LIS_constantsMod.o
read_ARMdata.o : LIS_fileIOMod.o
read_ARMdata.o : map_utils.o
read_ARMdata.o : LIS_logMod.o
read_ARMdata.o : LIS_coreMod.o
read_ARMdata.o : LIS_timeMgrMod.o
read_ARMdata.o : ARMdata_module.o
read_ARMdata.o : LIS_misc.h
