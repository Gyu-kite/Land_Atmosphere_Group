get_geos5fcst.o get_geos5fcst.d : get_geos5fcst.F90
get_geos5fcst.o : LDT_coreMod.o
get_geos5fcst.o : LDT_logMod.o
get_geos5fcst.o : LDT_metforcingMod.o
get_geos5fcst.o : LDT_constantsMod.o
get_geos5fcst.o : geos5fcst_forcingMod.o
get_geos5fcst.o : LDT_timeMgrMod.o
