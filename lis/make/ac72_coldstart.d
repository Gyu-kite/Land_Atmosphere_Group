ac72_coldstart.o ac72_coldstart.d : ac72_coldstart.F90
ac72_coldstart.o : LIS_coreMod.o
ac72_coldstart.o : LIS_timeMgrMod.o
ac72_coldstart.o : global.o
ac72_coldstart.o : LIS_logMod.o
ac72_coldstart.o : ac72_lsmMod.o
