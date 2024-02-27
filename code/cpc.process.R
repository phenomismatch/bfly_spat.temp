#Processing CPC climate data for use in analyses
#Elise Larsen, Georgetown U
#Run in R 4.1.2
#rain: CPC Global Unified Gauge-Based Analysis of Daily Precipitation data provided by the NOAA PSL, Boulder, Colorado, USA, from their website at https://psl.noaa.gov 
#yearly files: https://downloads.psl.noaa.gov/Datasets/cpc_global_precip/
#inputs yearly CPC data files (not included in repo)
#outputs: temp_mins_with_sd_09.RData, temp_maxs_with_sd_09.RData, precipSum.RData

#libraries

library(tidyverse)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(ggcorrplot)
library(RODBC)
library(lubridate)

#variables
eq_doys<-c(80,172,266,356)  #DOYs of solstices & equinoxes
save.output<-F #set to True to create output files from input files


##CPC Data Extraction

#Checking CPC dimensions
nc_data <- nc_open('data/cpc_data/tmax.2001.nc')
lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")
head(lon)

#function to extract temperature data by variable & year, with spatial hex grid
#change var to "tmin" to get minimum daily temps
cpc.temp.fx<- function(year, var="tmax") {
  filename=paste("data/cpc_data/",var,".",year,".nc", sep="")
  t1 <- brick(filename,varname=var,stopIfNotEqualSpaced=FALSE)
  e <- extent(200,350,0,70)
  t2<-crop(t1,e)
  
  #I needed to shift the x-axis to be on a -180 to 180 scale
  extent(t2)=c(xmin(t2)-360, xmax(t2)-360, ymin(t2), ymax(t2))

  hge <- rgdal::readOGR('data/maps/hex_grid_crop.shp', verbose = FALSE)
  hex <- spTransform(hge,projection(t2))
  r1.vals <- extract(t2, hex, fun=mean,na.rm=TRUE,df=TRUE,layer=1)
  r1.vals <- r1.vals %>%
    pivot_longer(
      cols = starts_with("X"),
      names_to = "Date",
      names_prefix = "X",
      values_to = paste("mean.",var,sep=""),
      values_drop_na = TRUE) %>%
    mutate(year=substr(Date,1,4),month=substr(Date,6,7),day=substr(Date,9,10), var=var)
  r2.vals <- extract(t2, hex, fun=sd,na.rm=TRUE,df=TRUE,layer=1)
  r2.vals <- r2.vals %>%
    pivot_longer(
      cols = starts_with("X"),
      names_to = "Date",
      names_prefix = "X",
      values_to = paste("sd."),
      values_drop_na = TRUE) %>%
    mutate(year=substr(Date,1,4),month=substr(Date,6,7),day=substr(Date,9,10), var=var)
  r.vals<-merge(r1.vals, r2.vals, by=c("ID","Date","year","month","day","var"))  
    
  r.vals$cell<-hex$cell[r.vals$ID]
  #r2.vals$cell<-hex$cell[r2.vals$ID]
  return(r.vals)
}

#function to extract precipitation data by variable & year, with spatial hex grid
cpc.pcp.fx<- function(year, var="precip") {
  filename=paste("data/cpc_data/",var,".",year,".nc", sep="")
  t1 <- brick(filename,varname=var,stopIfNotEqualSpaced=FALSE)
  e <- extent(200,350,0,70)
  t2<-crop(t1,e)
  
  #shift the x-axis to be on a -180 to 180 scale
  extent(t2)=c(xmin(t2)-360, xmax(t2)-360, ymin(t2), ymax(t2))
  
  hge <- rgdal::readOGR('data/maps/hex_grid_crop.shp', verbose = FALSE)
  hex <- spTransform(hge,projection(t2))
  r1.vals <- extract(t2, hex, fun=mean,na.rm=TRUE,df=TRUE,layer=1)
  r1.vals <- r1.vals %>%
    pivot_longer(
      cols = starts_with("X"),
      names_to = "Date",
      names_prefix = "X",
      values_to = paste("mean.",var,sep=""),
      values_drop_na = TRUE) %>%
    mutate(year=substr(Date,1,4),month=substr(Date,6,7),day=substr(Date,9,10), var=var)
  r2.vals <- extract(t2, hex, fun=sd,na.rm=TRUE,df=TRUE,layer=1)
  r2.vals <- r2.vals %>%
    pivot_longer(
      cols = starts_with("X"),
      names_to = "Date",
      names_prefix = "X",
      values_to = paste("sd.",var,sep=""),
      values_drop_na = TRUE) %>%
    mutate(year=substr(Date,1,4),month=substr(Date,6,7),day=substr(Date,9,10), var=var)
  r.vals<-bind_cols(r1.vals, r2.vals$sd.tmin)

  r.vals$cell<-hex$cell[r.vals$ID]

  return(r.vals)
}

## DATA EXTRACTION

##EXTRACT TMIN AND TMAX
tmin.list<-list()
tmax.list<-list()
for (y in 1:length(years)) {
  tmin.list[[y]]<-cpc.temp.fx(years[y],"tmin")
  tmax.list[[y]]<-cpc.temp.fx(years[y],"tmax")
}

tmin1<-as_tibble(bind_rows(tmin.list)) %>% dplyr::select(-var) %>% rename(sd.tmin=sd.)
tmax1<-as_tibble(bind_rows(tmax.list)) %>% dplyr::select(-var) %>% rename(sd.tmax=sd.)
if(save.output) {
  save(tmin1, file="temp_mins_with_sd_09.RData")
  save(tmax1, file="temp_maxs_with_sd_09.RData")
}

##EXTRACT PRECIPITATION
precip.l1<-list()
for (y in 1:length(years)) {
  precip.l1[[y]]<-cpc.pcp.fx(years[y],"precip")
}

##Extract precipitation data

pcp.data<-pcp.extract %>%
  filter(year>1998) %>%
  mutate(qtr=ceiling(as.numeric(month)/3)) %>%
  group_by(cell,year,qtr) %>% 
  summarize(qtr.precip=sum(mean.precip), nvals=tally())


##PROCESS CLIMATE DATA
#eq_doys<-c(80,172,266,356)
precip<-bind_rows(precip.l1) %>%
  mutate(date=as.Date(Date, format="%Y.%m.%d"), year=as.numeric(year), 
         month=as.numeric(month), day=as.numeric(day), doy=yday(date),
         yr.1=ifelse(month<10,year,year+1), 
         season=ifelse(doy<eq_doys[1],"winter.p",ifelse(doy<eq_doys[2],"spring.p",ifelse(doy<eq_doys[3],"summer.p","winter.p"))))
summary(precip)

precip.sum<-precip %>%
  group_by(cell,yr.1,season) %>%
  summarize(pcp.tot=sum(mean.precip, na.rm=T), pcp.sd=sd(mean.precip, na.rm=T)) %>%
  group_by(cell, season) %>%
  mutate(mean.pcp=median(pcp.tot, na.rm=T), mean.sd=median(pcp.sd, na.rm=T)) %>%
  ungroup() %>%
  group_by(cell,yr.1,season) %>%
  mutate(pcp.dev=(pcp.tot-mean.pcp), pcp.dev.std=(pcp.tot-mean.pcp)/mean.pcp, sd.dev=(pcp.sd-mean.sd)/mean.sd)
summary(precip.sum)  

if(save.output) {
save(precip.sum,file="data/derived/precipSum.RData")
}

cpc.data<-merge(tmin1, tmax1, by=intersect(names(tmin1),names(tmax1)))
nrow(cpc.data)
summary(cpc.data)
#


### CPC DERIVED METRICS 
#library(lubridate)
#load degree day function
source("code/src/degday1.R")

cpc.all<-cpc.data %>%
  dplyr::mutate(doy=yday(as.Date(paste(year,month, day, sep="-"),format="%Y-%m-%d")),
         season1=ifelse(doy>eq_doys[3],as.numeric(year)+1,as.numeric(year)),
         season2=ifelse(doy>eq_doys[4],4,ifelse(doy>eq_doys[3],3,ifelse(doy>eq_doys[2],2,ifelse(doy>eq_doys[1],1,4)))),
         warm=ifelse(mean.tmin>0,1,0), cold=ifelse(mean.tmax<=0,1,0), dailyGDD=0,  dailyGDDa=0,  dailyGDDb=0)

for(i in 1:nrow(cpc.all)) {
  cpc.all$dailyGDD[i]<-degreedays(cpc.all$mean.tmin[i], cpc.all$mean.tmax[i], 10,33)
}

cpc.metrics<-cpc.all %>%
  filter(season1 %in% (2000:2020)) %>%
  group_by(cell, season1, season2) %>%
  summarize(cold=sum(cold*floor(season2/2.5), na.rm=T),warm=sum(warm*floor(season2/2.5), na.rm=T),GDD=sum(dailyGDD, na.rm=T)) %>%
  group_by(cell, season1) %>%
  mutate(colddays=sum(cold, na.rm=T),warmdays=sum(warm, na.rm=T)) %>%
  filter(season2 %in% c(1,2))

if(save.output) {
  save(cpc.metrics, file="data/processed_cpc.RData")
}
