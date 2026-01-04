interp_geos5fcst.o interp_geos5fcst.d : interp_geos5fcst.F90
interp_geos5fcst.o : LDT_spatialDownscalingMod.o
interp_geos5fcst.o : LDT_coreMod.o
interp_geos5fcst.o : geos5fcst_forcingMod.o
