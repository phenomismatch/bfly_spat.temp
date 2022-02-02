#Hadley Project
#Butterfly Data import and summary
#Elise Larsen, started 2021-04-16

#libraries
library(tidyverse)
library(rinat)
library(ggplot2)
library(lubridate)
library(RODBC)
library(sp)
library(sf)
library(rgdal)
library(raster)


### Create grid

### Region of Interest
bounds <- c(25, -126, 52, -75)

#NABA data
con1 <- odbcConnect("NABA")
sqlTables(con1)

con2 <- odbcConnect("2019BFLY")
sqlTables(con2)

#######################  MIGRATION

##Species list
splist<-read_csv("species.list.csv")

sqlTables(con2)

spnames1<-sqlFetch(con2,"UniqueNamesResolution") %>%
  mutate(Species.std=paste(`NABA Genus`,`NABA Species`,sep=" ")) %>%
  filter(Species.std%in%splist$Species)
  
data1<-sqlFetch(con2,"1_Master_Obs_No dups")

data.filt<-merge(spnames1, data1, by.x="UMD_CODE", by.y="UMD-CODE")

table(data.filt$UMD_CODE)

data.locs<-sqlFetch(con2,"1_Master_Location_All no dups")

obs.data<-merge(data.filt, data.locs, by.x="GULocationID", by.y="GULocationID") %>%
  filter(between(GULongitude,bounds[2],bounds[4]) & between(GULatitude,bounds[1],bounds[3]))

obs.data<-obs.data %>%
  dplyr::select(UMD_CODE, Species.std, Program=Program.x, BasisOfRecord, ObsMonth, ObsDay, ObsYear, NumSeen, CountryClean, StateClean, GULatitude, GULongitude) %>%
  mutate(latid=floor(GULatitude*2)/2+0.25,lonid=floor(GULongitude*2)/2+0.25, ObsDate=as.Date(paste(ObsMonth,ObsDay,ObsYear,sep="/"),"%m/%d/%Y"), season=ifelse(ObsMonth<7,"spring","fall"))
  
obs.1<-na.omit(obs.data) %>%   mutate(DOY=yday(ObsDate),cellid=paste(lonid,latid,sep=".")) 


naba.summary<-obs.1 %>%
  group_by(Species.std, ObsYear, season, cellid) %>%
  tally()
naba.summary %>% filter(n>7) %>% group_by(Species.std) %>% tally()

naba.metrics<-obs.1 %>%
  group_by(lonid,latid,ObsYear,season,UMD_CODE) %>%
  add_tally() %>%
  filter(n>4) %>%
  summarize(medianDOY=median(DOY, na.rm=T), n.obs=mean(n, na.rm=T))

sp.doy<-naba.metrics %>%
  dplyr::select(lonid:medianDOY) %>%
  pivot_wider(names_from="UMD_CODE", values_from="medianDOY")
names(sp.doy)[5:10]<-c("VANCAR.median","DANPLE.median","AGRAVAN.median","EURPRO.median","CHICAT.median","ASCMON.median")

sp.n<-naba.metrics %>%
  dplyr::select(lonid:UMD_CODE,n.obs) %>%
  pivot_wider(names_from="UMD_CODE", values_from="n.obs")
names(sp.n)[5:10]<-c("VANCAR.n","DANPLE.n","AGRAVAN.n","EURPRO.n","CHICAT.n","ASCMON.n")

metrics<-merge(sp.doy,sp.n,  by=intersect(names(sp.doy),names(sp.n)))

write.csv(metrics, file="mig_metrics.5d_alldata.csv")


########### Mean of first n(n=3 or 5) observations
spnames<-c("VANCAR","DANPLE","AGRAVAN","EURPRO","CHICAT","ASCMON")


nmin<-5
naba.metrics<-obs.1 %>%
  group_by(lonid,latid,ObsYear,season,UMD_CODE) %>%
  dplyr::mutate(v=1, n1=cumsum(v), ntot=max(n1)) %>%
  #add_tally() %>%
  filter(ntot>4, n1<=nmin) %>%
  summarize(medianDOY=median(DOY, na.rm=T), meanDOY=mean(DOY, na.rm=T), n.obs=mean(ntot, na.rm=T))


################################### OLD



  


#fetch observations table, filter to surveys with at least 10 species in the applicable years
naba.obs<-sqlFetch(con1,"2_NABA_sightings_flatfile") %>%
  mutate(spname=paste(`Gen/Tribe/Family`,Species,sep=" ")) %>%
  filter(spname %in% splist$Species, between(UMD_Latitude, bounds[1],bounds[3]), between(UMD_Longitude, bounds[2],bounds[4])) %>%
  mutate(latid=floor(UMD_Latitude),lonid=floor(UMD_Longitude),cellid=paste(lonid,latid,sep="."))

#
naba.obs<-sqlFetch(con1,"NABA_OBS_FlatFile") %>%
  mutate(spname=paste(`Gen/Tribe/Fam`,Species_Epithet,sep=" ")) %>%
  filter(spname %in% splist$Species, between(Latitude, bounds[1],bounds[3]), between(Longitude, bounds[2],bounds[4])) %>%
  mutate(latid=floor(Latitude*2)/2+0.25,lonid=floor(Longitude*2)/2+0.25,cellid=paste(lonid,latid,sep="."))


naba.summary<-naba.obs %>%
  group_by(spname, Year, cellid) %>%
  tally()
naba.summary %>% filter(n>7) %>% group_by(spname) %>% tally()

#naba.sp<-st_as_sf(naba.obs, coords = c("UMD_Longitude", "UMD_Latitude"), crs = 4326) %>%
#  st_transform(crs = 3857)

table(naba.obs$spname, naba.obs$Year)

summary(naba.obs)

#INaturalist data
inat.extract<-list()

for(i in 1:nrow(splist)) {
  inat.extract[[i]]<-get_inat_obs(query = splist$Species[i], bounds=bounds) %>%
    mutate(obsdate=as.Date(datetime), year=year(obsdate),lonid=floor(longitude), latid=floor(latitude), cellid=paste(lonid, latid, sep="."))
  
}

save(inat.extract, naba.obs, file="rawdata202108.RData")


allinat<-as.tibble(bind_rows(inat.extract)) %>%
  mutate(obsdate=as.Date(datetime), year=year(obsdate),lonid=floor(longitude*2)/2+0.5, latid=floor(latitude*2)/2+0.5, cellid=paste(lonid, latid, sep=".")) %>%
  filter(quality_grade=="research") %>%
  dplyr::select(cellid, year, latid, lonid, obsdate, spname=scientific_name)
naba1<-naba.obs %>% dplyr::select(cellid, latid, lonid, year=Year, obsdate=SurveyDate, spname)

obs.data<-bind_rows(allinat, naba1)
obs.data.summary<-obs.data %>%
  group_by(spname, year, cellid) %>%
  tally()
obs.data.summary %>% filter(n>7) %>% group_by(spname) %>% tally()

inatnames<-sort(unique(allinat$spname))
matched<-c("Ascia monuste","Ascia monuste","Calycopis isobean",NA,"Danaus plexippus",
           "Agraulis vanillae","Agraulis vanillae",NA,"Eurema proterpia","Vanessa cardui")

obs.data<-obs.data %>% mutate(name=spname)
for(i in 1:nrow(obs.data)) {
  obs.data$name[i]<-matched[which(inatnames==obs.data$spname[i])]
}

first.pass<-obs.data %>%
  mutate(doy=yday(obsdate), season=ifelse(doy>183,"fall","spring")) %>%
  group_by(name, latid, lonid, year, season) %>%
  add_tally() %>%
  filter(n>4) %>%
  summarize(median=median(doy))

pheno1<-first.pass %>%
  mutate(cell.latitude=latid-0.25, cell.longitude=lonid-0.25)
  
  
pheno1<-pheno1[,c(3:8)]
write.csv(pheno1, file="bfly_metrics_draft.csv")

#amonuste<-get_inat_obs(query = "Ascia monuste", bounds=bounds)
#amonuste<-amonuste %>%
#  mutate(obsdate=as.Date(datetime), lonid=floor(longitude*2)/2, latid=floor(latitude*2)/2, cellid=paste(lonid, latid, sep="."))
#hist(year(amonuste$obsdate))
#table(amonuste$cellid, year(amonuste$obsdate))
