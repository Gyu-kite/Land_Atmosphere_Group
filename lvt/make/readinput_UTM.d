readinput_UTM.o readinput_UTM.d : readinput_UTM.F90
readinput_UTM.o : LVT_domainMod.o
readinput_UTM.o : LVT_logMod.o
readinput_UTM.o : LVT_coreMod.o
readinput_UTM.o : map_utils.o
readinput_UTM.o : LVT_misc.h
