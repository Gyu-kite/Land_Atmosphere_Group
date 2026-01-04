read_netcdf_4d.o read_netcdf_4d.d : read_netcdf_4d.F90
read_netcdf_4d.o : LDT_coreMod.o
read_netcdf_4d.o : LDT_gridmappingMod.o
read_netcdf_4d.o : LDT_misc.h
read_netcdf_4d.o : CLM45_parmsMod.o
