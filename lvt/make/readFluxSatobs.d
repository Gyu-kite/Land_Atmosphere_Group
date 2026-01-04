readFluxSatobs.o readFluxSatobs.d : readFluxSatobs.F90
readFluxSatobs.o : FluxSat_obsMod.o
readFluxSatobs.o : LVT_timeMgrMod.o
readFluxSatobs.o : LVT_histDataMod.o
readFluxSatobs.o : LVT_constantsMod.o
readFluxSatobs.o : LVT_coreMod.o
readFluxSatobs.o : LVT_misc.h
readFluxSatobs.o : LVT_logMod.o
