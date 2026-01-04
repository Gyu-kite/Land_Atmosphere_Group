get_mrms_grib.o get_mrms_grib.d : get_mrms_grib.F90
get_mrms_grib.o : LIS_coreMod.o
get_mrms_grib.o : mrms_grib_forcingMod.o
get_mrms_grib.o : LIS_logMod.o
get_mrms_grib.o : LIS_constantsMod.o
get_mrms_grib.o : LIS_timeMgrMod.o
