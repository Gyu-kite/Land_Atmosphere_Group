lnd_atmMod.o lnd_atmMod.d : lnd_atmMod.F90
lnd_atmMod.o : LIS_precisionMod.o
lnd_atmMod.o : LIS_misc.h
lnd_atmMod.o : clm2_varcon.o
lnd_atmMod.o : clm2_varsur.o
lnd_atmMod.o : clm2_lsmMod.o
lnd_atmMod.o : clm2_varpar.o
