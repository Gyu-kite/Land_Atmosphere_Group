interp_geos5fcst.o interp_geos5fcst.d : interp_geos5fcst.F90
interp_geos5fcst.o : LIS_spatialDownscalingMod.o
interp_geos5fcst.o : geos5fcst_forcingMod.o
interp_geos5fcst.o : LIS_coreMod.o
