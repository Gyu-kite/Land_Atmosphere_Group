finalize_stg2.o finalize_stg2.d : finalize_stg2.F90
finalize_stg2.o : stg2_forcingMod.o
finalize_stg2.o : LDT_coreMod.o
finalize_stg2.o : LDT_misc.h
