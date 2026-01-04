read_HYMAP_basin.o read_HYMAP_basin.d : read_HYMAP_basin.F90
read_HYMAP_basin.o : LDT_logMod.o
read_HYMAP_basin.o : LDT_coreMod.o
read_HYMAP_basin.o : HYMAP_parmsMod.o
read_HYMAP_basin.o : LDT_fileIOMod.o
