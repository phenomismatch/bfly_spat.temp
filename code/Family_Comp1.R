#NABA BIS butterfly family composition

#libraries
library(tidyverse)
library(RODBC)
library(lubridate)
library(readxl)
library(stringr)


##Import NABA Data
#Database connection
con1 <- odbcConnect("NABA")
sqlTables(con1)

naba.bis<-sqlFetch(con1,"1_BISSightings all records") 
summary(naba.bis)
naba.obs<-sqlFetch(con1,"NABA_OBS_FlatFile") 
summary(naba.obs)
naba.sp<-sqlFetch(con1,"0_NABA BFLY NAMES") 
summary(naba.sp)
naba.lh<-sqlFetch(con1,"Life history table for Allen") 
summary(naba.lh)
write.csv(naba.lh, file="fam.comp1.csv")

0_NABA BFLY NAMES

temp_PhenologyMacro1D_summary
NABA_OBS_FlatFile 
Life history table for Allen
1_BISSightings all records
2_NABA_sightings_flatfile
