get_geos5fcst.o get_geos5fcst.d : get_geos5fcst.F90
get_geos5fcst.o : LIS_timeMgrMod.o
get_geos5fcst.o : LIS_metforcingMod.o
get_geos5fcst.o : geos5fcst_forcingMod.o
get_geos5fcst.o : LIS_coreMod.o
get_geos5fcst.o : LIS_logMod.o
get_geos5fcst.o : LIS_constantsMod.o
