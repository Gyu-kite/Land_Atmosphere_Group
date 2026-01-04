read_WRFAK_elev.o read_WRFAK_elev.d : read_WRFAK_elev.F90
read_WRFAK_elev.o : LDT_misc.h
read_WRFAK_elev.o : LDT_metforcingMod.o
read_WRFAK_elev.o : WRF_AKdom_forcingMod.o
read_WRFAK_elev.o : LDT_fileIOMod.o
read_WRFAK_elev.o : LDT_coreMod.o
read_WRFAK_elev.o : LDT_logMod.o
