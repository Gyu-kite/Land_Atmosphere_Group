interp_petusgs.o interp_petusgs.d : interp_petusgs.F90
interp_petusgs.o : LIS_coreMod.o
interp_petusgs.o : LIS_logMod.o
interp_petusgs.o : petusgs_forcingMod.o
