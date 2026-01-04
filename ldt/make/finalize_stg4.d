finalize_stg4.o finalize_stg4.d : finalize_stg4.F90
finalize_stg4.o : stg4_forcingMod.o
finalize_stg4.o : LDT_misc.h
finalize_stg4.o : LDT_coreMod.o
