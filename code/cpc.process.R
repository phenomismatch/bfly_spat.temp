#rain: CPC Global Unified Gauge-Based Analysis of Daily Precipitation data provided by the NOAA PSL, Boulder, Colorado, USA, from their website at https://psl.noaa.gov 
#yearly files: https://downloads.psl.noaa.gov/Datasets/cpc_global_precip/

#libraries

library(tidyverse)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(ggcorrplot)
library(RODBC)
#
eq_doys<-c(80,172,266,356)

                   
#nc_rain_data <- nc_open('C:/Users/eal109/Downloads/precip.1999.nc') #precip at 0.5 degree scale
#nc_rain_data <- nc_open('C:/Users/eal109/Downloads/precip.mon.total.v7.nc') #precip at 0.5 degree scale
#lon <- ncvar_get(nc_rain_data, "lon")
#lat <- ncvar_get(nc_rain_data, "lat", verbose = F)
#t <- ncvar_get(nc_rain_data, "time")
#precip.array <- ncvar_get(nc_data, "precip") # store the data in a 3-dimensional array
#dim(precip.array)


##CPC TEMPERATURE
nc_data <- nc_open('C:/Users/eal109/Documents/Macrosystems/PredictCaterpillarPheno/cpc_data/tmax.2001.nc')
lon <- ncvar_get(nc_data, "lon")

lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")
head(lon)

tmax.array <- ncvar_get(nc_data, "tmax") # store the data in a 3-dimensional array
dim(tmax.array)

cpc.temp.fx<- function(year, var="tmax") {
  filename=paste("C:/Users/eal109/Documents/Macrosystems/PredictCaterpillarPheno/cpc_data/",var,".",year,".nc", sep="")
  t1 <- brick(filename,varname=var,stopIfNotEqualSpaced=FALSE)
  e <- extent(200,350,0,70)
  t2<-crop(t1,e)
  
  #I needed to shift the x-axis to be on a -180 to 180 scale
  extent(t2)=c(xmin(t2)-360, xmax(t2)-360, ymin(t2), ymax(t2))

  hge <- rgdal::readOGR('C:/Users/eal109/Documents/Bird_Phenology/Data/hex_grid_crop/hex_grid_crop.shp', verbose = FALSE)
  hex <- spTransform(hge,projection(t2))
  r.vals <- extract(t2, hex, fun=mean,na.rm=TRUE,df=TRUE,layer=1)
  r.vals <- r.vals %>%
    pivot_longer(
      cols = starts_with("X"),
      names_to = "Date",
      names_prefix = "X",
      values_to = paste("mean.",var,sep=""),
      values_drop_na = TRUE) %>%
    mutate(year=substr(Date,1,4),month=substr(Date,6,7),day=substr(Date,9,10), var=var)
  
  r.vals$cell<-hex$cell[r.vals$ID]
  return(r.vals)
}

cpc.pcp.fx<- function(year, var="precip") {
  #filename=('C:/Users/eal109/Downloads/precip.mon.total.v7.nc') #precip at 0.5 degree scale
  filename=paste("C:/Users/eal109/Documents/Macrosystems/PredictCaterpillarPheno/cpc_data/",var,".",year,".nc", sep="")
  t1 <- brick(filename,varname=var,stopIfNotEqualSpaced=FALSE)
  e <- extent(200,350,0,70)
  t2<-crop(t1,e)
  
  #shift the x-axis to be on a -180 to 180 scale
  extent(t2)=c(xmin(t2)-360, xmax(t2)-360, ymin(t2), ymax(t2))
  
  hge <- rgdal::readOGR('C:/Users/eal109/Documents/Bird_Phenology/Data/hex_grid_crop/hex_grid_crop.shp', verbose = FALSE)
  hex <- spTransform(hge,projection(t2))
  r.vals <- extract(t2, hex, fun=mean,na.rm=TRUE,df=TRUE,layer=1)
  r.vals <- r.vals %>%
    pivot_longer(
      cols = starts_with("X"),
      names_to = "Date",
      names_prefix = "X",
      values_to = paste("mean.",var,sep=""),
      values_drop_na = TRUE) %>%
    mutate(year=substr(Date,1,4),month=substr(Date,6,7),day=substr(Date,9,10), var=var)
  
  r.vals$cell<-hex$cell[r.vals$ID]
  return(r.vals)
}



#pcp.extract<-cpc.pcp.fx(var="precip")
pcp.data<-pcp.extract %>%
  filter(year>1998) %>%
  mutate(qtr=ceiling(as.numeric(month)/3)) %>%
  group_by(cell,year,qtr) %>% 
  summarize(qtr.precip=sum(mean.precip), nvals=tally())

pcp.data<-pcp.data %>%
  mutate(season=ifelse(qtr<4,as.numeric(year),as.numeric(year)-1))

years<-c(1999:2020)
tmin.list<-list()
tmax.list<-list()
precip.list<-list()
for (y in 1:length(years)) {
  tmin.list[[y]]<-cpc.temp.fx(years[y],"tmin")
  tmax.list[[y]]<-cpc.temp.fx(years[y],"tmax")
}

precip.l1<-list()
for (y in 1:length(years)) {
  precip.l1[[y]]<-cpc.pcp.fx(years[y],"precip")
}
library(lubridate)

load("data/derived/precip.RData")
eq_doys<-c(80,172,266,356)
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
save(precip.sum,file="data/derived/precipSum.RData")
 
tmin1<-as_tibble(bind_rows(tmin.list)) %>% dplyr::select(-var)
tmax1<-as_tibble(bind_rows(tmax.list)) %>% dplyr::select(-var)
cpc.data<-merge(tmin1, tmax1, by=intersect(names(tmin1),names(tmax1)))
nrow(cpc.data)
summary(cpc.data)
#


### CPC DERIVED METRICS 
library(lubridate)
#load degree day function
source("C:/Users/eal109/Downloads/Git/Git/LabTraining/01GitRStudio/src/degday1.R")

cpc.all<-cpc.data %>%
  dplyr::mutate(doy=yday(as.Date(paste(year,month, day, sep="-"),format="%Y-%m-%d")),
         season1=ifelse(doy>eq_doys[3],as.numeric(year)+1,as.numeric(year)),
         season2=ifelse(doy>eq_doys[4],4,ifelse(doy>eq_doys[3],3,ifelse(doy>eq_doys[2],2,ifelse(doy>eq_doys[1],1,4)))),
         warm=ifelse(mean.tmin>0,1,0), cold=ifelse(mean.tmax<=0,1,0), dailyGDD=0)

for(i in 1:nrow(cpc.all)) {
  cpc.all$dailyGDD[i]<-degreedays(cpc.all$mean.tmin[i], cpc.all$mean.tmax[i], 10,33)
}
save(cpc.all, file="cpc_data/processed_cpc.RData")

cpc.metrics<-cpc.all %>%
  filter(season1 %in% (2000:2020)) %>%
  group_by(cell, season1, season2) %>%
  summarize(cold=sum(cold*floor(season2/2.5), na.rm=T),warm=sum(warm*floor(season2/2.5), na.rm=T),GDD=sum(dailyGDD, na.rm=T)) %>%
  group_by(cell, season1) %>%
  mutate(colddays=sum(cold, na.rm=T),warmdays=sum(warm, na.rm=T)) %>%
  filter(season2 %in% c(1,2))

sp1<-cpc.metrics %>%
  filter(season2==1) %>%
  mutate(spring.gdd = GDD) %>%
  dplyr::select(year=season1, cell, colddays, warmdays, spring.gdd)

sp2<-cpc.metrics %>%
  filter(season2==2) %>%
  mutate(summer.gdd = GDD) %>%
  dplyr::select(year=season1, cell, colddays, warmdays, summer.gdd)

therm<-merge(sp1, sp2, by=intersect(names(sp1), names(sp2)))
summary(therm)
write.csv(therm, file="cpc_data/therm_metrics.csv")

