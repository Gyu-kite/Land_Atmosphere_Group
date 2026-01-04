read_geos5fcst.o read_geos5fcst.d : read_geos5fcst.F90
read_geos5fcst.o : LDT_coreMod.o
read_geos5fcst.o : LDT_logMod.o
read_geos5fcst.o : LDT_misc.h
read_geos5fcst.o : geos5fcst_forcingMod.o
read_geos5fcst.o : LDT_metforcingMod.o
