#Elise Larsen, Georgetown U
##Accepted 2023-12
##Examining variation due to GDD thresholds and within-hex variation

#libraries
library(tidyverse)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(corrplot)
library(ggcorrplot)
library(RODBC)
library(lubridate)

#days to delineate seasons in "season2" field
eq_doys<-c(80,172,266,356)
#between eq_doy[1] and eq_doy[2] is N. hemisphere spring = 1
#between eq_doy[2] and eq_doy[3] is N. hemisphere summer = 2
#between eq_doy[3] and eq_doy[4] is N. hemisphere autumn = 3
#below eq_doy[1] or above eq_doy[4] is N. hemisphere winter = 4

##CPC TEMPERATURE
nc_data <- nc_open('data/cpc_data/tmax.2001.nc')
lon <- ncvar_get(nc_data, "lon")

lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")
head(lon)

#This function captures both mean and sd of daily temperature metrics
cpc.temp.fx2<- function(year, var="tmax") {
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
  return(r.vals)
}

years<-c(2010:2019)
tmin.list<-list()
tmax.list<-list()

for (y in 1:length(years)) {
  tmin.list[[y]]<-cpc.temp.fx2(years[y],"tmin")
  tmax.list[[y]]<-cpc.temp.fx2(years[y],"tmax")
}


tmin1<-as_tibble(bind_rows(tmin.list)) %>% dplyr::select(-var) %>% rename(sd.tmin=sd.)
tmax1<-as_tibble(bind_rows(tmax.list)) %>% dplyr::select(-var) %>% rename(sd.tmax=sd.)
save(tmin1, file="temp_mins_with_sd_10.RData")
save(tmax1, file="temp_maxs_with_sd_10.RData")

cpc.data<-merge(tmin1, tmax1, by=intersect(names(tmin1),names(tmax1)))
nrow(cpc.data)
summary(cpc.data)
#


### CPC DERIVED METRICS 
library(lubridate)
#load degree day function
source("code/src/degday1.R")

cpc.all<-cpc.data %>%
  dplyr::mutate(doy=yday(as.Date(paste(year,month, day, sep="-"),format="%Y-%m-%d")),
                season1=ifelse(doy>eq_doys[3],as.numeric(year)+1,as.numeric(year)),
                season2=ifelse(doy>eq_doys[4],4,ifelse(doy>eq_doys[3],3,ifelse(doy>eq_doys[2],2,ifelse(doy>eq_doys[1],1,4)))),
                warm=ifelse(mean.tmin>0,1,0), cold=ifelse(mean.tmax<=0,1,0), dailyGDD=0,  dailyGDDa=0,  dailyGDDb=0)


##Calculating daily GDD values with different thresholds, different tmin & tmax
for(i in 1:nrow(cpc.all)) {
  #calculate GDD as used in the main analysis
  cpc.all$dailyGDD[i]<-degreedays(cpc.all$mean.tmin[i], cpc.all$mean.tmax[i], 10,33)
  #calculate GDD using lower 95% CI values for tmin and tmax
  cpc.all$dailyGDDsdlo[i]<-degreedays(cpc.all$mean.tmin[i]-(1.96*cpc.all$sd.tmin[i]), cpc.all$mean.tmax[i]-(1.96*cpc.all$sd.tmax[i]),10,33)
  #calculate GDD using higher 95% CI values for tmin and tmax
  cpc.all$dailyGDDsdhi[i]<-degreedays(cpc.all$mean.tmin[i]+(1.96*cpc.all$sd.tmin[i]), cpc.all$mean.tmax[i]+(1.96*cpc.all$sd.tmax[i]),10,33)
  #calculate GDD using higher 95% CI values for tmin and lower 95% CI values for tmax
  cpc.all$dailyGDDsdmin[i]<-degreedays(cpc.all$mean.tmin[i]+(1.96*cpc.all$sd.tmin[i]), cpc.all$mean.tmax[i]-(1.96*cpc.all$sd.tmax[i]),10,33)
  #calculate GDD using lower 95% CI values for tmin and higher 95% CI values for tmax
  cpc.all$dailyGDDsdmax[i]<-degreedays(cpc.all$mean.tmin[i]-(1.96*cpc.all$sd.tmin[i]), cpc.all$mean.tmax[i]+(1.96*cpc.all$sd.tmax[i]),10,33)
  #calculate GDD using lower temperature thresholds by 3 degrees
  cpc.all$dailyGDDa[i]<-degreedays(cpc.all$mean.tmin[i], cpc.all$mean.tmax[i], 7,30)
  #calculate GDD using higher temperature thresholds by 3 degrees
  cpc.all$dailyGDDb[i]<-degreedays(cpc.all$mean.tmin[i], cpc.all$mean.tmax[i], 13,36)
}


#filter to cells in study area
load("data/spatial.domain.RData")
#filters from 141 cells to 54 of 68 in study region
#and calculates seasonal summary GDD values by hex & year
## Spring GDD (from DOY 80 to DOY 172: approximating spring equinox to summer solstice)
cpc.metrics<-filter(cpc.all, cell %in% STUDYCELLS) %>%
  filter(season1 %in% (2000:2019) & season2 %in% c(1,2)) %>%
  group_by(cell, season1, season2) %>%
  summarize(GDD=sum(dailyGDD, na.rm=T),GDD0730=sum(dailyGDDa, na.rm=T),GDD1336=sum(dailyGDDb, na.rm=T),
            GDDsd.lo=sum(dailyGDDsdlo, na.rm=T),GDDsd.hi=sum(dailyGDDsdhi, na.rm=T),
            GDDsd.min=sum(dailyGDDsdmin, na.rm=T),GDDsd.max=sum(dailyGDDsdmax, na.rm=T)) %>%
  rename(year=season1)
 
#save(cpc.metrics, file="data/processed_cpc_newthresholds_revtest.RData")

#Calculate annual deviations

#create recent baseline
env_baseline<-cpc.metrics %>%
  mutate(base.year=ifelse(year>2015,1,NA)) %>%
  mutate(gdd=GDD*base.year,gdd07=GDD0730*base.year,gdd13=GDD1336*base.year,
         gddlow=GDDsd.lo*base.year, gddhigh=GDDsd.hi*base.year,
         gddmin=GDDsd.min*base.year, gddmax=GDDsd.max*base.year) %>%
  group_by(cell,season2) %>%
  mutate(base.gdd=mean(gdd, na.rm=T), base.gdd07=mean(gdd07, na.rm=T),base.gdd13=mean(gdd13, na.rm=T),
         base.gddlow=mean(gddlow, na.rm=T), base.gddhigh=mean(gddhigh, na.rm=T),
         base.gddmin=mean(gddmin, na.rm=T), base.gddmax=mean(gddmax, na.rm=T)       ) 

env_deviations<-env_baseline %>%
  mutate(gdd.dev=gdd-base.gdd, gdd.dev07=gdd07-base.gdd07, gdd.dev13=gdd13-base.gdd13,
         gdd.low.dev=gddlow-base.gddlow, gdd.high.dev=gddhigh-base.gddhigh,
         gdd.min.dev=gddmin-base.gddmin, gdd.max.dev=gddmax-base.gddmax  ) %>%
  dplyr::select(cell, year, season2, gdd.dev:gdd.max.dev)


env1<-filter(env_deviations, !is.na(gdd.dev))
corrplot(cor(env1[,4:10]))
cor(env1[,4:10])

save(env1, file="data/derived/gdd.variations.RData")

#Function to extract metrics from pairwise regressions

result.table<-c("1",0,1,1,0)

extract_cor<-function(lmsum) {
  result<-c(lmsum$coefficients[1,1],lmsum$coefficients[2,1],lmsum$adj.r.squared[1], lmsum$sigma)
  return(result)
  }
t1<-env1 %>%
  dplyr::select(gdd.dev07:gdd.max.dev) %>%
  map_dfr(function(x) extract_cor(summary(lm(env1$gdd.dev~x))))

t2<-as.data.frame(t(as.matrix(t1)))


blah.df<-as.data.frame("int"=t1[1,],"coef"=t1[2,], "adj.r2"=t1[3,])


unlist(t1)

t2<-data.frame(matrix(unlist(t1),byrow=T))


save(cpc.all, file="data/processed_cpc_rev10.RData")

## SUMMER CLIMATE

#load("data/processed_cpc_rev.RData")
#load("data/processed_cpc_newthresholds_rev.RData")

for(i in 1:nrow(cpc.all)) {
  #cpc.all$dailyGDD[i]<-degreedays(cpc.all$mean.tmin[i], cpc.all$mean.tmax[i], 10,33)
  cpc.all$dailyGDDlow[i]<-degreedays(cpc.all$mean.tmin[i]-cpc.all$sd.tmin[i], cpc.all$mean.tmax[i]-cpc.all$sd.tmax[i],10,33)
  cpc.all$dailyGDDhigh[i]<-degreedays(cpc.all$mean.tmin[i]+cpc.all$sd.tmin[i], cpc.all$mean.tmax[i]+cpc.all$sd.tmax[i],10,33)
  cpc.all$dailyGDDc[i]<-degreedays(cpc.all$mean.tmin[i], cpc.all$mean.tmax[i], 5,30)
}

cpc.metrics2<-cpc.all %>%
  filter(season1 %in% (2000:2020)) %>%
  group_by(cell, season1, season2) %>%
  summarize(cold=sum(cold*floor(season2/2.5), na.rm=T),warm=sum(warm*floor(season2/2.5), na.rm=T),GDD=sum(dailyGDD, na.rm=T),
            GDDlow=sum(dailyGDDlow, na.rm=T),GDDhigh=sum(dailyGDDhigh, na.rm=T),GDD5=sum(dailyGDDc, na.rm=T)) %>%
  group_by(cell, season1) %>%
  mutate(colddays=sum(cold, na.rm=T),warmdays=sum(warm, na.rm=T)) %>%
  filter(season2 %in% c(1,2))
save(cpc.metrics2, file="data/derived/processed_cpc_newthresholds.RData")
#load("data/derived/processed_cpc_newthresholds.RData")



summary(lm(GDDlow~GDD, data=cpc.metrics2))
summary(lm(GDDhigh~GDD, data=cpc.metrics2))

ggplot(data=cpc.metrics2, aes(x=GDD, y=GDDlow)) + geom_point(color="blue") + geom_point(aes(x=GDD, y=GDDhigh),color="red")  + geom_point(aes(x=GDD, y=GDD5), color="green") + labs(y="GDD") + 
  theme_minimal()


sp1<-cpc.metrics2 %>%
  filter(season2==1) %>%
  mutate(spring.gdd = GDD, spring.low = GDDlow, spring.high = GDDhigh, spring.GDD5=GDD5) %>%
  dplyr::select(year=season1, cell, colddays, warmdays, spring.gdd, spring.low, spring.high, spring.GDD5)

sp2<-cpc.metrics2 %>%
  filter(season2==2) %>%
  mutate(summer.gdd = GDD, summer.low = GDDlow, summer.high = GDDhigh, summer.GDD5=GDD5) %>%
  dplyr::select(year=season1, cell, colddays, warmdays, summer.gdd, summer.low, summer.high, summer.GDD5)

therm<-merge(sp1, sp2, by=intersect(names(sp1), names(sp2)))
summary(therm)
write.csv(therm, file="data/derived/therm_metrics_rev.csv")

#plot correlations
corrplot(cor(sp1[,5:8]))
cor(sp1[,5:8])
corrplot(cor(sp1[,5:8]))
cor(sp2[,5:8])

clim.data<-therm %>%
  filter(cell %in% STUDYCELLS)


##
env_baseline<-clim.data %>%
  filter(year<2020) %>%
  mutate(base.year=ifelse(year>2015,1,NA)) %>%
  mutate(spring=spring.gdd*base.year,summer=summer.gdd*base.year,
         spring5=spring.GDD5*base.year,summer5=summer.GDD5*base.year) %>%
  group_by(cell) %>%
  mutate(base.spring=mean(spring, na.rm=T), base.summer=mean(summer, na.rm=T),
         base.spring5=mean(spring5, na.rm=T), base.summer5=mean(summer5, na.rm=T)) %>%
  mutate(spring.gdd.dev=spring.gdd-base.spring, summer.gdd.dev=summer.gdd-base.summer,
         spring.gdd5.dev=spring.GDD5-base.spring5, summer.gdd5.dev=summer.GDD5-base.summer5,
         spring.gddlow.dev=spring.low-base.spring, summer.GDDlow.dev=summer.low-base.summer,
         spring.gddhigh.dev=spring.high-base.spring, summer.gddhigh.dev=summer.high-base.summer) %>%
  select(year:cell,spring.gdd.dev:summer.gddhigh.dev)



