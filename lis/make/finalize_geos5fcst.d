finalize_geos5fcst.o finalize_geos5fcst.d : finalize_geos5fcst.F90
finalize_geos5fcst.o : geos5fcst_forcingMod.o
finalize_geos5fcst.o : LIS_coreMod.o
finalize_geos5fcst.o : LIS_misc.h
