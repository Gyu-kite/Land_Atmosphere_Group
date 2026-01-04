write_lookup_table.o write_lookup_table.d : write_lookup_table.F90
write_lookup_table.o : LDT_coreMod.o
write_lookup_table.o : LDT_logMod.o
write_lookup_table.o : SMOSNRTNNL2sm_obsMod.o
write_lookup_table.o : LDT_constantsMod.o
write_lookup_table.o : LDT_NetCDF_inc.h
write_lookup_table.o : LDT_misc.h
