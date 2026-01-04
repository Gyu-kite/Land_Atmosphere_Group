finalize_ecmwf.o finalize_ecmwf.d : finalize_ecmwf.F90
finalize_ecmwf.o : LIS_coreMod.o
finalize_ecmwf.o : LIS_misc.h
finalize_ecmwf.o : ecmwf_forcingMod.o
