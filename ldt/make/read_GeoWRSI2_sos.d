read_GeoWRSI2_sos.o read_GeoWRSI2_sos.d : read_GeoWRSI2_sos.F90
read_GeoWRSI2_sos.o : LDT_coreMod.o
read_GeoWRSI2_sos.o : LDT_logMod.o
read_GeoWRSI2_sos.o : GeoWRSI_parmsMod.o
read_GeoWRSI2_sos.o : map_utils.o
read_GeoWRSI2_sos.o : geowrsi2_arraymgmt_module.o
read_GeoWRSI2_sos.o : fbil_module.o
