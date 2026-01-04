finalize_WRFout.o finalize_WRFout.d : finalize_WRFout.F90
finalize_WRFout.o : LIS_coreMod.o
finalize_WRFout.o : WRFout_forcingMod.o
