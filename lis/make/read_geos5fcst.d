read_geos5fcst.o read_geos5fcst.d : read_geos5fcst.F90
read_geos5fcst.o : LIS_constantsMod.o
read_geos5fcst.o : LIS_misc.h
read_geos5fcst.o : LIS_coreMod.o
read_geos5fcst.o : LIS_logMod.o
read_geos5fcst.o : geos5fcst_forcingMod.o
read_geos5fcst.o : LIS_metforcingMod.o
