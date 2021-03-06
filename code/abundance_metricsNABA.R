#NABA Abundance Metrics
## Resident butterfly species groups defined by overwinter stage (OWS)
## Abundances = natural log of summed butterfly per hour abundances within OWS group
## NABA counts from "seasonal surveys" (raw data not publically available)
## Species OWS traits compiled by GU Ries Lab
#E Larsen, Georgetown U, Updated 2021-05


#libraries
library(tidyverse)
library(RODBC)
library(lubridate)
library(readxl)
library(stringr)

##Set parameters
minSR<-10 #minimum recorded species richness (across all species) to include a survey
maxBPH<-100 #maximum butterflies per hour value to include; records above this are thrown out
month.filt<-c(6:8) #months in which to include surveys
year.filt<-c(2000:2017) #years in which to include surveys
minYears<-2 #minimum number of years with survey data for a count circle to be included

##Import NABA Data

#Import NABA count circle locations, mapped to hex cell
naba.circles<-read_csv("data/naba/naba_circles.csv") 

#Database connection
con1 <- odbcConnect("NABAcircle2018")
sqlTables(con1)

#fetch observations table, filter to surveys with at least 10 species in the applicable years
naba.obs<-sqlFetch(con1,"NFJ_Observations to 2018")  %>% 
  group_by(CountID, SurveyID) %>% add_tally(name="RSR") %>% 
  filter(RSR>=minSR) %>% #, ObsYear%in%year.filt) %>%
  dplyr::select(ID:SightID,FromDate:NumSeen)

#fetch surveys table
#filter by year and month, Surveys meeting minimum reported species richness, Circles surveyed in a min # years
naba.surveys<-sqlFetch(con1,"NFJ_Surveys to 2018") %>% 
  dplyr::select(CountID, SurveyID, SurveyDate, ObsYear, ObsMonth,Lat,Lng,Party_Hours) %>%
  mutate(doy=yday(SurveyDate)) %>%
  filter(SurveyID %in% naba.obs$SurveyID, ObsMonth %in% month.filt) %>% 
  group_by(CountID) %>%  mutate(nYears=length(unique(ObsYear))) %>%
  filter(CountID %in% naba.circles$CountID, nYears>=minYears)

#Merge surveys and circle locations to tag by hex cell
survey.data<-merge(naba.surveys, naba.circles, by.x=c("CountID","Lat","Lng"), by.y=c("CountID","Lat","Lng")) 

##For Casey
casey.naba<-survey.data %>%
  dplyr::select(CountID, SurveyID, Lat, Lng, ObsYear) %>%
  dplyr::group_by(Lat, Lng, ObsYear) %>%
  tally()
  #dplyr::summarize(minyr=min(ObsYear, na.rm=T), maxyr=max(ObsYear, na.rm=T), nyr=n()) 
casey.naba

write.csv(casey.naba, file="nabaAbund_locYear.csv")
##Import Species Trait Data & Extract OWS group + Species 
trait.data<-read_xlsx("data/naba/SpeciesTraitsLR.xlsx") %>% 
  mutate(ScientificName=paste(`NABA Genus`,ifelse(is.na(`NABA species`),"",`NABA species`),sep=" ")) %>%
  dplyr::select(ID, ScientificName, group=`Simpleton grouping code`) %>%
  mutate(ScientificName=str_trim(ScientificName, side="both")) %>%
  filter(group %in% c("RE","RL","RP")) %>%
  group_by(ScientificName, group) %>% tally()

#Merge observation and trait data, match sci names using naba.names table
naba.names<-read_csv(file="data/naba/naba_names.csv") %>% dplyr::select(ScientificName, SpeciesID)
naba.obs<-merge(na.omit(naba.obs), naba.names,by=c("ScientificName"))  %>%
  dplyr::select(CountID:ObsYear,Lat:NumSeen,ScientificName,SpeciesID) %>% 
  filter(SurveyID %in% survey.data$SurveyID) 

obs.data<-merge(naba.obs, trait.data, by.x=c("SpeciesID"), by.y=c("ScientificName"))

##Calculate abundance metrics by survey for OWS species groups
#Calculate butterflies per party hour (bph) for each species recorded
#Sum bph values within OWS species groups, log-transform, & round values to thousandths 
#report observed species richness within groups
naba.abundance<-merge(obs.data, survey.data, by=intersect(names(obs.data), names(survey.data))) %>%
  mutate(bph=NumSeen/Party_Hours, doy=yday(SurveyDate), sp1=1) 

(naba.abundance %>% group_by(SpeciesID) %>% add_tally(name = 'nrec') %>% filter(bph>maxBPH) %>% group_by(SpeciesID) %>% summarize(n100=n(), perc=round(n100/nrec,4)) %>% group_by(SpeciesID) %>% summarize(n100=mean(n100, na.rm=T), perc=mean(perc, na.rm=T)))

naba.abundance<-merge(obs.data, survey.data, by=intersect(names(obs.data), names(survey.data))) %>%
mutate(bph=NumSeen/Party_Hours, doy=yday(SurveyDate), sp1=1) %>%
group_by(cell, ObsYear, CountID,SurveyID, SurveyDate) %>% mutate(SR=sum(sp1)) %>%
filter(bph<=maxBPH, SR>=minSR) %>%  ##, SR>=minSR
group_by(cell, Lat, Lng, CountID, SurveyID, ObsYear, ObsMonth, doy, group ) %>%
summarize(abund.bph=round(sum(bph),3), log.abund=round(log(sum(bph)),3), SR=sum(sp1))


#Write abundance metrics data table to csv
write.csv(naba.abundance, file="data/derived/naba_OWS_abundances.csv")

naba.abundance %>% filter(abund.bph>100)
summary(naba.abundance)

naba.abundance<-naba.abundance %>%
  dplyr::select(cell,latitude=Lat, longitude=Lng, CountID, year=ObsYear, doy, group, abund.bph, n.sp=SR)
write.csv(naba.abundance, file="NABA_grouped_abundance.csv")

naba1<-naba.abundance %>% filter(Lat<30,abund.bph<100, group=="RL") %>% 
  group_by(CountID,Lat,Lng) %>% tally() %>%
  group_by(CountID) %>% tally()

#summary(table(naba1$CountID))

