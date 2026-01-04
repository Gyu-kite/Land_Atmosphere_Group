get_ncep_filenames.o get_ncep_filenames.d : get_ncep_filenames.F90
get_ncep_filenames.o : LIS_logMod.o
get_ncep_filenames.o : LIS_coreMod.o
get_ncep_filenames.o : LIS_forecastMod.o
