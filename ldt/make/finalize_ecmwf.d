finalize_ecmwf.o finalize_ecmwf.d : finalize_ecmwf.F90
finalize_ecmwf.o : ecmwf_forcingMod.o
finalize_ecmwf.o : LDT_misc.h
finalize_ecmwf.o : LDT_coreMod.o
