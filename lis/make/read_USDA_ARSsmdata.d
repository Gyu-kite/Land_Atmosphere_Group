read_USDA_ARSsmdata.o read_USDA_ARSsmdata.d : read_USDA_ARSsmdata.F90
read_USDA_ARSsmdata.o : LIS_fileIOMod.o
read_USDA_ARSsmdata.o : map_utils.o
read_USDA_ARSsmdata.o : USDA_ARSsm_obsMod.o
read_USDA_ARSsmdata.o : LIS_timeMgrMod.o
read_USDA_ARSsmdata.o : LIS_coreMod.o
read_USDA_ARSsmdata.o : LIS_logMod.o
read_USDA_ARSsmdata.o : LIS_mpiMod.o
read_USDA_ARSsmdata.o : LIS_constantsMod.o
