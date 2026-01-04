finalize_AWRAL.o finalize_AWRAL.d : finalize_AWRAL.F90
finalize_AWRAL.o : LIS_misc.h
finalize_AWRAL.o : AWRAL_forcingMod.o
finalize_AWRAL.o : LIS_coreMod.o
