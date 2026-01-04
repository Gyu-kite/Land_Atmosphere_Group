finalize_climatology.o finalize_climatology.d : finalize_climatology.F90
finalize_climatology.o : LIS_coreMod.o
finalize_climatology.o : climatology_VariablesMod.o
finalize_climatology.o : LIS_misc.h
