readcrd_mrms_grib.o readcrd_mrms_grib.d : readcrd_mrms_grib.F90
readcrd_mrms_grib.o : LIS_logMod.o
readcrd_mrms_grib.o : LIS_coreMod.o
readcrd_mrms_grib.o : mrms_grib_forcingMod.o
