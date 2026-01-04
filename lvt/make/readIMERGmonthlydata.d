readIMERGmonthlydata.o readIMERGmonthlydata.d : readIMERGmonthlydata.F90
readIMERGmonthlydata.o : LVT_coreMod.o
readIMERGmonthlydata.o : IMERG_monthly_dataMod.o
readIMERGmonthlydata.o : LVT_constantsMod.o
readIMERGmonthlydata.o : LVT_timeMgrMod.o
readIMERGmonthlydata.o : LVT_misc.h
readIMERGmonthlydata.o : LVT_histDataMod.o
readIMERGmonthlydata.o : LVT_logMod.o
