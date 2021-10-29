#Environmental variables for phenology & abundance modeling
#E Larsen, Georgetown U, Updated 2021-05


##Sources
#GDD, frost-free days and frozen days calcs use CPC data:
#Cite: CPC Global Temperature data provided by the NOAA/OAR/ESRL PSL, Boulder, Colorado, USA, from their Web site at https://psl.noaa.gov/data/gridded/data.cpc.globaltemp.html

#libraries
library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library(sp)

#Set Parameters
eq_doys<-c(80,172,266,356) #DOYs for equinox & solstice dates
#cells in study area
load("data/spatial.domain.RData")


## SUMMER CLIMATE
#Calculate summary seasonal GDD values for each hex cell & year:
## Spring GDD (from DOY 80 to DOY 172: approximating spring equinox to summer solstice)

#CPC data has spring gdd, summer gdd, winter cold days (tmax<=0), and winter warm days (tmin>0)
clim.data<-read_csv("data/derived/therm_metrics.csv") %>%
  filter(cell %in% STUDYCELLS)

## GREENUP
###Timing of half-max greenup (forest pixels)

#Greenup: halfmax DOY forest pixels
green.for<-readRDS("data/envir/MidGreenup-2020-08-06-forest.rds") %>%
  dplyr::select(year, cell, gr_mn_for=gr_mn, gr_pfor=gr_pcell, cell_lat, cell_lng) %>%
  filter(cell %in% STUDYCELLS, !is.na(gr_mn_for))


#Greenup: half-max DOY all pixels
green.all<-readRDS("data/envir/MidGreenup-2020-08-06-all.rds") %>% 
  dplyr::select(year, cell, gr_mn_all=gr_mn, gr_pall=gr_pcell, cell_lat, cell_lng) %>%
  filter(cell %in% STUDYCELLS, !is.na(gr_mn_all))


#Greenup: calculate half-max DOY of open pixels
greendev<-merge(green.all,green.for, by=intersect(names(green.all),names(green.for)), all.x=T, all.y=T) %>%
  mutate(p_open=gr_pall-gr_pfor) %>%
  mutate(gr_mn_open=ifelse(p_open>0,((gr_mn_all - gr_mn_for*gr_pfor)/(gr_pall-gr_pfor)),NA), gr_mn_lag=(gr_mn_open-gr_mn_for))


env_var<-merge(greendev, clim.data, by=intersect(names(greendev), names(clim.data)), all.x=T, all.y=T)

ggplot(data=filter(greendev, year>2009), aes(x=gr_mn_for, y=gr_mn_lag, color=gr_pfor)) + geom_point() + 
  labs(x="Greenup in forest pixels", y="Lag of open pixels", color="percent forest") + 
  theme_minimal()

write.csv(env_var, file="data/derived/environmental.vars.csv")



#### PCA

## filter data to what is needed and PCA 
load("data/derived/pheno.RData")
abund<-read_csv("data/derived/naba_OWS_abundances.csv")%>%
  dplyr::select(cell, year=ObsYear) %>%
  summarize(cellyr=paste(cell,year,sep=".")) 
cells.ab<-sort(unique(abund$cellyr))
cells.ph<-sort(unique(paste(pheno.quant$cell,pheno.quant$year,sep=".")))

env_var<-env_var %>% mutate(cellyr=paste(cell,year,sep=".")) %>%
  filter(cellyr %in% c(cells.ab, cells.ph))

corr <- round(cor(na.omit(env_var[,c(7,11,15)])), 1)
ggcorrplot(corr, method = "circle")


templm<-lm(year~gr_mn_for+gr_mn_lag+spring.gdd+summer.gdd, data=env_var)
library(car)
vif(templm)


pca_res <- prcomp(na.omit(env_var[,c(7,15)]), scale. = TRUE)
summary(pca_res)
pca_res
pca<-data.frame(cell=env_var$cell[(!is.na(env_var$spring.gdd) & !is.na(env_var$gr_mn_for))], year=env_var$year[(!is.na(env_var$spring.gdd) & !is.na(env_var$gr_mn_for))], pc1=(pca_res$x[,1]),pc2=-pca_res$x[,2] )


env2<-merge(env_var,pca, by=intersect(names(env_var),names(pca)),all.x=T)
corr <- round(cor(env2[,c(11,13,14,16,18,19)]), 1)
ggcorrplot(corr, method = "circle")

write.csv(env2, file="data/derived/envir2.csv")
####

### PLOT
ggplot(data=env2, aes(x=year, y=spring.gdd)) + geom_point()

