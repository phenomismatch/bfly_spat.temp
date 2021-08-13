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

## SUMMER CLIMATE
#Calculate summary seasonal GDD values for each hex cell & year:
## Spring GDD (from DOY 80 to DOY 172: approximating spring equinox to summer solstice)
#Daily mean growing degree day values by hex cell
cpc.gdd<-na.omit(read_csv("data/envir/cpcGDD_hex_2000-2020.csv")) %>%
  dplyr::select(cell=hex_cell, date, year, doy, meanGDD) %>%
  filter(doy>=eq_doys[1],doy<=eq_doys[3]) %>%
  mutate(season=ifelse(doy<eq_doys[2],"spring.gdd","summer.gdd")) %>%
  group_by(cell, year, season) %>%
  summarize(gdd=sum(meanGDD, na.rm=T)) %>%
  pivot_wider(names_from = season, values_from = gdd) %>%
  group_by(cell) %>% mutate(gdd.comb=spring.gdd+summer.gdd, spring.mean=mean(spring.gdd, na.rm=T),spring.dev=(spring.gdd-spring.mean),
                            summer.mean=mean(summer.gdd, na.rm=T),summer.dev=(summer.gdd-summer.mean)) 


## WINTER CLIMATE




cpc.winter<-read.csv("data/envir/cpc_winter.csv" ) %>%
  dplyr::select(cell=hex_cell, year=season, FRD, FFD) %>%
  group_by(cell) %>%
  filter(year %in% c(2001:2020)) %>%
  mutate(FFDhm=mean(FFD, na.rm=T),FRDhm=mean(FRD, na.rm=T)) %>%
  group_by(cell, year) %>%
  summarize(FRD=mean(FRD, na.rm=T), FFD=mean(FFD, na.rm=T),FFD.dev=ifelse(is.na(FFD),NA,(FFD-FFDhm)),FRD.dev=ifelse(is.na(FRD),NA,(FRD-FRDhm))) %>%
  dplyr::select(cell, year, FFD, FRD, FFD.dev, FRD.dev)


clim.data<-merge(cpc.gdd, cpc.winter, by=c("cell","year"), all.x=T, all.y=T)

write.csv(clim.data, file="data/envir/cpc.comb.csv")
## GREENUP
###Timing of half-max greenup (forest pixels)

#Greenup: halfmax DOY forest pixels
green.for<-readRDS("data/envir/MidGreenup-2020-08-06-forest.rds") %>%
  dplyr::select(year, cell, gr_mn_for=gr_mn, gr_pfor=gr_pcell, cell_lat, cell_lng)

#Greenup: half-max DOY all pixels
green.all<-readRDS("data/envir/MidGreenup-2020-08-06-all.rds") %>% 
  dplyr::select(year, cell, gr_mn_all=gr_mn, gr_pall=gr_pcell, cell_lat, cell_lng)

#Greenup: calculate half-max DOY of open pixels
greendev<-merge(green.all,green.for, by=intersect(names(green.all),names(green.for)), all.x=T, all.y=T) %>%
  mutate(p_open=gr_pall-gr_pfor) %>%
  mutate(gr_mn_open=ifelse(p_open>0,((gr_mn_all*gr_pall - gr_mn_for*gr_pfor)/(gr_pall-gr_pfor)),NA), gr_mn_lag=(gr_mn_open-gr_mn_for))

env_var<-merge(greendev, clim.data, by=intersect(names(greendev), names(clim.data)), all.x=T, all.y=T)

ggplot(data=filter(greendev, year>2009), aes(x=gr_mn_for, y=gr_mn_lag, color=gr_pfor)) + geom_point() + 
  labs(x="Greenup in forest pixels", y="Lag of open pixels", color="percent forest") + 
  theme_minimal()

write.csv(env_var, file="data/derived/environmental.vars.csv")

## filter data to what is needed and PCA 
abund<-read_csv("data/derived/naba_OWS_abundances.csv") %>%
  dplyr::select(cell, year=ObsYear) %>%
  summarize(cellyr=paste(cell,year,sep=".")) 

tempabund<-abund %>% dplyr::select(cellyr) %>% group_by(cellyr) %>% tally()
pheno.datafile<-"data/derived/adult_bfly_metrics.csv"
pheno.orig<-read_csv(pheno.datafile)
hexyrs<-pheno.orig %>% dplyr::select(year, cell=HEXcell) %>%
  mutate(cellyr=paste(cell,year,sep=".")) %>% 
  group_by(year, cell, cellyr) %>% tally()

env_var<-env_var %>% mutate(cellyr=paste(cell,year,sep=".")) %>%
  filter(cellyr %in% c(tempabund$cellyr,hexyrs$cellyr))

corr <- round(cor(na.omit(env_var[,c(7,10,11,14,16,18,21:22)])), 1)
ggcorrplot(corr, method = "circle")


pca_res <- prcomp(na.omit(env_var[,c(7,10,14)]), scale. = TRUE)
summary(pca_res)
pca_res
pca<-data.frame(cell=env_var$cell[!is.na(env_var$gr_mn_for)], year=env_var$year[!is.na(env_var$gr_mn_for)],  warmearly=(pca_res$x[,1]),warmlateopen=-pca_res$x[,2] )


env2<-merge(env_var[,-23],pca, by=intersect(names(env_var),names(pca)),all.x=T)
corr <- round(cor(env2[,c(23,24,16,18,21,22)]), 1)
ggcorrplot(corr, method = "circle")

write.csv(env2, file="data/derived/envir.csv")
