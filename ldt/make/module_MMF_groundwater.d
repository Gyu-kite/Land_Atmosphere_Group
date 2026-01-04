module_MMF_groundwater.o module_MMF_groundwater.d : module_MMF_groundwater.F90
module_MMF_groundwater.o : LDT_logMod.o
module_MMF_groundwater.o : LDT_coreMod.o
module_MMF_groundwater.o : LDT_paramMaskCheckMod.o
module_MMF_groundwater.o : LDT_constantsMod.o
module_MMF_groundwater.o : LDT_gridmappingMod.o
module_MMF_groundwater.o : map_utils.o
module_MMF_groundwater.o : LDT_paramDataMod.o
module_MMF_groundwater.o : LDT_misc.h
