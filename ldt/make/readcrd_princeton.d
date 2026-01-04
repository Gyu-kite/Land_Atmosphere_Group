readcrd_princeton.o readcrd_princeton.d : readcrd_princeton.F90
readcrd_princeton.o : princeton_forcingMod.o
readcrd_princeton.o : LDT_coreMod.o
readcrd_princeton.o : LDT_logMod.o
