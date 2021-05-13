#NABA Abundance Metrics
## Resident butterfly species groups defined by overwinter stage (OWS)
## Abundances = natural log of summed butterfly per hour abundances within OWS group
## NABA counts from "seasonal surveys" 
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
month.filt<-c(6,7) #months in which to include surveys
year.filt<-c(2000:2018) #years in which to include surveys
minYears<-2 #minimum number of years with survey data for a count circle to be included

##Import NABA Data

#Import NABA count circle locations, mapped to hex cell
naba.circles<-read_csv("adult_data/naba_circles.csv") 

#Database connection
con1 <- odbcConnect("NABAcircle2018")

#fetch observations table, filter to surveys with at least 10 species in the applicable years
naba.obs<-sqlFetch(con1,"NFJ_Observations to 2018")  %>% 
  group_by(CountID, SurveyID) %>% add_tally(name="RSR") %>% 
  filter(RSR>=minSR, ObsYear%in%year.filt) %>%
  select(ID:SightID,FromDate:NumSeen)

#fetch surveys table
#filter by year and month, Surveys meeting minimum reported species richness, Circles surveyed in a min # years
naba.surveys<-sqlFetch(con1,"NFJ_Surveys to 2018") %>% 
  select(CountID, SurveyID, SurveyDate, ObsYear, ObsMonth,Lat,Lng,Party_Hours) %>%
  mutate(doy=yday(SurveyDate)) %>%
  filter(SurveyID %in% naba.obs$SurveyID, ObsMonth %in% month.filt) %>% 
  group_by(CountID) %>%  mutate(nYears=length(unique(ObsYear))) %>%
  filter(CountID %in% naba.circles$CountID, nYears>=minYears)

#Merge surveys and circle locations to tag by hex cell
survey.data<-merge(naba.surveys, naba.circles, by.x=c("CountID","Lat","Lng"), by.y=c("CountID","Lat","Lng")) 


##Import Species Trait Data & Extract OWS group + Species 
trait.data<-read_xlsx("adult_data/SpeciesTraitsLR.xlsx") %>% 
  mutate(ScientificName=paste(`NABA Genus`,ifelse(is.na(`NABA species`),"",`NABA species`),sep=" ")) %>%
  select(ID, ScientificName, group=`Simpleton grouping code`) %>%
  mutate(ScientificName=str_trim(ScientificName, side="both")) %>%
  filter(group %in% c("RE","RL","RP")) %>%
  group_by(ScientificName, group) %>% tally()

#Merge observation and trait data, match sci names using naba.names table
naba.names<-read_csv(file="naba_names.csv") %>% select(ScientificName, SpeciesID)
naba.obs<-merge(na.omit(naba.obs), naba.names,by=c("ScientificName"))  %>%
  select(CountID:ObsYear,Lat:NumSeen,ScientificName,SpeciesID) %>% 
  filter(SurveyID %in% survey.data$SurveyID) 

obs.data<-merge(naba.obs, trait.data, by.x=c("SpeciesID"), by.y=c("ScientificName"))

naba.abundance<-merge(obs.data, survey.data, by=intersect(names(obs.data), names(survey.data))) %>%
  mutate(bph=NumSeen/Party_Hours, doy=yday(SurveyDate), sp1=1) %>%
  group_by(SurveyID) %>% mutate(SR=sum(sp1)) %>%
  filter(bph<=maxBPH, SR>=minSR) %>%
  group_by(cell, Lat, Lng, CountID, SurveyID, ObsYear, ObsMonth, doy, group) %>%
  summarize(abund.bph=sum(bph), log.abund=log(abund.bph), SR=sum(sp1))

write.csv(naba.abundance, file="data/derived_data/naba_OWS_abundances.csv")




























#NABA Abundance Metrics
#E Larsen, Georgetown U, 2021

#libraries
library(tidyverse)
library(RODBC)
library(sp)
library(lubridate)
library(readxl)
library(ggplot2)
library(mrfDepth)
library(stringr)

#Import NABA Data
con1 <- odbcConnect("NABAcircle2018")
#tbls <- sqlTables(con1)
#tbls$TABLE_NAME
naba.circle.obs<-sqlFetch(con1,"NFJ_Observations to 2018")

#Filter to surveys with at least 10 species
naba.obs<-naba.circle.obs %>% group_by(CountID, SurveyID) %>% add_tally(name="RSR") %>% filter(RSR>9)

naba.surveys<-sqlFetch(con1,"NFJ_Surveys to 2018") %>% 
  select(CountID, SurveyID, SurveyDate, ObsYear, ObsMonth,Lat,Lng,Party_Hours) %>%
  mutate( doy=yday(SurveyDate)-min(yday(SurveyDate)+1))

naba.circles<-read_csv("data/naba/naba_circles.csv")

#Filter to surveys in June & July, 2000-2017, with location and abundance data
surveys.filtered<-naba.surveys %>%
  filter(CountID %in% naba.circles$CountID, SurveyID %in% naba.obs$SurveyID, ObsMonth %in% c(6,7), ObsYear %in% 2000:2017)
  
#Survey circles with at least 2 years of data
naba.surveys<-merge(surveys.filtered, naba.circles, by.x=c("CountID","Lat","Lng"), by.y=c("CountID","Lat","Lng")) %>%
  group_by(CountID) %>%  mutate(nYears=length(unique(ObsYear))) %>%
  filter( nYears>1)



traits<-read_xlsx("data/naba/SpeciesTraitsLR.xlsx")
#table(traits$`Simpleton grouping code`)

traits<-traits %>% 
  mutate(ScientificName=paste(`NABA Genus`,ifelse(is.na(`NABA species`),"",`NABA species`),sep=" ")) %>%
  select(ID, ScientificName, group=`Simpleton grouping code`) %>%
  mutate(ScientificName=str_trim(ScientificName, side="both")) %>%
  group_by(ScientificName, group) %>% tally()


naba.sp<-naba.obs %>% 
  group_by(ScientificName) %>%
  tally() %>%
  mutate(sciname=str_trim(ScientificName,side="both"), l1=str_count(sciname, " ")) %>%
  mutate(SpeciesID=sciname)

naba.sp<-na.omit(naba.sp)
for(rowi in 1:nrow(naba.sp)) {
  if(naba.sp$l1[rowi]>1) {
    naba.sp$SpeciesID[rowi]<-substr(naba.sp$sciname[rowi],1, (gregexpr(" ", naba.sp$sciname[rowi])[[1]][2]-1))
  }
}
naba.names<-sort(unique(naba.sp$SpeciesID))
unmatched<-naba.names[is.na(match(naba.names, traits$ScientificName))]


naba.obs<-na.omit(naba.obs) %>% select(CountID:ObsYear,Lat:NumSeen) %>% filter(SurveyID %in% naba.surveys$SurveyID) 

naba.obs<-merge(naba.obs, naba.sp[,c(1,5)], by.x=c("ScientificName"), by.y=c("ScientificName"))

obs.data<-merge(naba.obs, traits, by.x=c("SpeciesID"), by.y=c("ScientificName"))


naba1<-merge(obs.data, naba.surveys, by=intersect(names(obs.data), names(naba.surveys))) %>%
  mutate(bph=NumSeen/Party_Hours, doy=yday(SurveyDate), sp1=1) %>%
  group_by(SurveyID) %>% mutate(SR=sum(sp1)) %>%
  filter(bph<=100, SR>9) 

naba.abundance<-naba1 %>%
  filter(group %in% c("RE","RL","RP")) %>%
  group_by(cell, Lat, CountID, SurveyID, ObsYear, ObsMonth, doy, group) %>%
  summarize(abund.bph=sum(bph), log.abund=log(abund.bph), SR=sum(sp1))

save(naba.abundance, file="naba_OWS_abundances.RData")
write.csv(naba.abundance, file="naba_OWS_abundances.csv")

naba.tot<-naba1 %>%
  group_by(cell, Lat, CountID, SurveyID,ObsYear, ObsMonth, doy) %>%
  summarize(abund.bph=sum(bph), log.abund=log(abund.bph), groupSR=sum(sp1))

ggplot(data=naba.tot, aes(group=ObsYear, y=log.abund)) + geom_boxplot() 


test<-lm(log.abund~Lat+ObsYear, data=naba.tot)
test<-lm(log.abund~Lat+ObsYear+group, data=naba.abundance)


##DEMOG VARIABLES
naba.groups.PY<-naba.groups %>% 
  select(cell, Lat, CountID, SurveyID, ObsYear,ObsMonth, doy, group, prev.ab=logall) %>%
  mutate(ObsYear=ObsYear+1) 
  

naba.data<-merge(naba.groups, naba.groups.PY[,c(1,3,5,8,9)], by=c("cell","CountID","ObsYear","group"), all.x=T)



###############################
##  ENVIRONMENTAL VARIABLES
load('pheno202105.RData')


##Combine tables
ab.data<-merge(naba.groups, naba.groups.PY, by=intersect(names(naba.groups),names(naba.groups.PY)))



ab.data2<-merge(ab.data, dataset2, by.x=c("cell","ObsYear"), by.y=c("HEXcell","year")) %>%
  mutate(sg=`Simpleton grouping code`, RL.log.abund=logall) %>% filter(sg=="RL")


library(lme4)

#model1<-lmer(logall~1+onset.dev+prev.ab+PC1+FR.dev+1|CountID, data=ab.data2)

length(unique(ab.data2$cell))
table(ab.data2$ObsYear)

ab.lm.full<-lm(RL.log.abund~1+onset.dev+duration.dev+prev.ab+FR.dev+as.factor(ObsMonth)+warmearly, data=ab.data2)
summary(ab.lm.full)
AIC(ab.lm.full)
pvals<-summary(ab.lm.full)$coefficients[c(2:nrow(summary(ab.lm.full)$coefficients)),4]
temp<-which(pvals)
plot.full<-plot_model(ab.lm.full) 
plot.full + theme_minimal() + geom_hline(yintercept=0) + 
  annotate(geom="text",y=c(rep(-.8,6)), x=c(6,5,4,3,2,1), label=ifelse(as.character(round(pvals,3))=="0","<0.01",as.character(round(pvals,3))))


ggplot(data=ab.data2, aes(x=ObsYear, y=RL.log.abund, color=as.factor(cell))) + geom_point()

lm2<-lm(RL.log.abund~1+as.factor(cell), data=ab.data2)
lm3<-lm(RL.log.abund~1+as.factor(ObsYear), data=ab.data2)



ab.lm.2<-lm(logall~1+onset.dev+duration.dev+prev.ab+ObsMonth+warmearly, data=ab.data2)
summary(ab.lm.2)
AIC(ab.lm.2)

plot_model(ab.lm.2)


ab.lmm.full<-lmer(logall~1+onset.dev+duration.dev+prev.ab+as.factor(ObsMonth)+warmearly+1|CountID, data=ab.data2)
summary(ab.lmm.full)
ab.lmm.x<-lmer(logall~1+prev.ab+warmearly+1|CountID, data=ab.data2)
summary(ab.lmm.x)



corrz2 <- round(cor(env2[,c("cell_lat","warmearly","coldearly","greenDZ","gdd_z2")]), 1)
ggcorrplot(corrz2, method = "circle",  lab = TRUE)






nabas<-naba1 %>% group_by(cell, ObsYear, CountID) %>% summarize(nct=1) %>% 
  group_by(cell, ObsYear) %>% summarize(cts=sum(nct, na.rm=T))
naba3<-nabas %>% group_by(cell) %>% tally()
summary(nabas)
length(unique(nabas$cell))

naba33<-nabas %>% filter(cts>2)

ggplot(data=naba1, aes(x=doy, y=logall, color=as.factor(ObsYear))) + geom_point() + theme(legend.position ="none") + 
  facet_wrap(~cell)

ggplot(data=filter(naba1, Lat<35), aes(x=doy, y=logall, color=as.factor(ObsYear))) + geom_point() + theme(legend.position ="none") + 
  labs(x="DOY",y="Log BPH", title="<35N") + theme_classic()

ggplot(data=filter(naba1, Lat>42), aes(x=doy, y=logall, color=as.factor(ObsYear))) + geom_point() + theme(legend.position ="none") + 
  labs(x="DOY",y="Log BPH", title=">42N") + theme_classic()

ggplot(data=filter(naba1, Lat<35), aes(x=ObsMonth, y=logall, group=ObsMonth)) + geom_boxplot() + theme(legend.position ="none") + 
  labs(x="Month",y="Log BPH", title="<35N") + theme_classic()

ggplot(data=filter(naba1, Lat>42), aes(x=ObsMonth, y=logall, group=ObsMonth)) + geom_boxplot() + theme(legend.position ="none") + 
  labs(x="MOnth",y="Log BPH", title=">42N") + theme_classic()


ggplot(data=filter(naba1, Lat>42, cell==534), aes(x=doy, y=logall, group=as.factor(cell), color=as.factor(cell))) + geom_point() + theme(legend.position ="none") + 
  labs(x="DOY",y="Log BPH", title=">42N") + theme_classic() + facet_wrap(~ObsYear)


naba2<-naba1 %>% 
  mutate(latrange=ifelse(Lat>42,"high",ifelse(Lat<35,"low","mid"))) %>%
   group_by(latrange, ObsYear, ObsMonth) %>% tally()

ggplot(data=filter(naba2, latrange=="low"), aes(x=ObsMonth, y=n, group=ObsMonth)) + geom_boxplot() + theme(legend.position ="none") + 
  labs(x="Month",y="N Surveys", title="<35N") + theme_classic()

ggplot(data=filter(naba2, latrange=="high"), aes(x=ObsMonth, y=n, group=ObsMonth)) + geom_boxplot() + theme(legend.position ="none") + 
  labs(x="MOnth",y="N Surveys", title=">42N") + theme_classic()



ggplot(data=filter(naba1, Lat<35), aes(x=doy, y=logall, color=as.factor(ObsYear))) + geom_point() + theme(legend.position ="none")

library(mgcv)
summary(gam(logall~doy + doy^2 +as.factor(ObsYear)+as.factor(cell), data=naba1))

naba.stats<-naba.obs %>% select(CountID, SurveyID, Lat, Lng, ObsYear, FromDate, ScientificName, NumSeen, `Simpleton grouping code`) %>%
  mutate(obsdoy=yday(FromDate), obsmonth=month(FromDate))
n1<-naba.obs %>% mutate(ObsDoy=yday(FromDate),ObsMonth=month(FromDate)) %>% 
  group_by(Lat,Lng,ObsYear,ObsMonth,ObsDoy,`Simpleton grouping code`) %>% #tally()
  summarize(tot=sum(NumSeen, na.rm=T)) %>% filter(between(ObsMonth,5,9))
ggplot(data=n1, aes(x=factor(ObsMonth), y=log(tot),fill=factor(`Simpleton grouping code`))) + 
  geom_boxplot() + labs(x="Month",y="Log Abundance",fill="Species group") + theme_classic()

hist((naba.stats$obsmonth))
ggplot(data=naba.stats, aes(x=as.factor(obsmonth), y=log(NumSeen))) + geom_boxplot()





naba.cs<-naba.circle.surveys %>%
  filter(CountID %in% naba.cl$CountID, ObsYear>1999)



