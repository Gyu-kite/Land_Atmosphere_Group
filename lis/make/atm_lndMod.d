atm_lndMod.o atm_lndMod.d : atm_lndMod.F90
atm_lndMod.o : lnd_atmMod.o
atm_lndMod.o : LIS_misc.h
atm_lndMod.o : LIS_precisionMod.o
atm_lndMod.o : LIS_timeMgrMod.o
atm_lndMod.o : clm2_shr_const_mod.o
