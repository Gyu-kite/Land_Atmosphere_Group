finalize_princeton.o finalize_princeton.d : finalize_princeton.F90
finalize_princeton.o : LIS_coreMod.o
finalize_princeton.o : princeton_forcingMod.o
finalize_princeton.o : LIS_misc.h
