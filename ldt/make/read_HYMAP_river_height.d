read_HYMAP_river_height.o read_HYMAP_river_height.d : read_HYMAP_river_height.F90
read_HYMAP_river_height.o : HYMAP_parmsMod.o
read_HYMAP_river_height.o : LDT_fileIOMod.o
read_HYMAP_river_height.o : LDT_logMod.o
read_HYMAP_river_height.o : LDT_coreMod.o
