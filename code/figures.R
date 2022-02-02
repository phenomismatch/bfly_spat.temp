##Larsen et al.
##Butterfly phenology & abundance by overwinter stage
##Submitted 2022-02
##Code by Elise Larsen, Georgetown University
##File contents: R code to produce manuscript figures (Figure 1 & some supplements)

#load libraries
library(ggeffects)
library(gridExtra)
library(lme4)
library(sf)
library(viridis)
library(tidyverse)
theme_set(theme_classic(base_size = 15))


# Import Data -------------------------------------------------------------

#file names: spatial domain & map features
spatial.domain<-"data/spatial.domain.RData"
hex.map.file<-"data/maps/hex_grid_crop.shp"
land.map.file<-"data/maps/ne_50m_admin_1_states_provinces_lakes.shp"

#file names: phenology & abundance metrics
pheno.datafile<-"data/derived/adult_bfly_phenometrics_noCountCircles_withFull2020Data.csv"
pheno.dev.file<-"data/derived/newphenodev.RData"
abundance.file<-"data/derived/naba_OWS_abundances.csv"
#other parameters
temporal.domain<-c(2001:2017)
#Map figure axis limits
f1lims.x<-c(-10700000,-6679169)
f1lims.y<-c(3500000,6500000)
ows.colors<-viridis_pal()(8)[c(1,4,7)]

#spatial domain (hex grid cells)
load(spatial.domain)  ## STUDYCELLS vector

pt1 = st_point(c(-100,24))
pt2 = st_point(c(-60,50))

bounds<-st_sfc(pt1, pt2) # %>% st_transform(crs ="+init=epsg:3857")
sf_points <- 
  bounds %>% 
  st_as_sf(coords = c('a', 'b'), crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs ="+init=epsg:3857")
box1<-st_bbox(sf_points)

## POLYGONS FOR MAP FIGURES
#load polygons of study hexes, reproject, and filter to study cells (#3857 is projected pseudo-mercator, 4326 is lat/long)
hex_sf <- st_transform(read_sf(hex.map.file), 3857) %>%   filter(cell %in% STUDYCELLS)
summary(hex_sf)

nam_sf <- st_transform(read_sf(land.map.file) , 3857) %>%
  filter(sr_adm0_a3 %in% c("CAN", "USA"), iso_3166_2 != "US-AK", iso_3166_2 != "US-HI") 


#import & filter phenology data
pheno.quant<-read_csv(pheno.datafile) %>% rename(cell=HEXcell) %>%
  filter(cell %in% STUDYCELLS, year %in% temporal.domain, between(q50,152,243), !is.na(q5))
#calculate # years of data per cell
pheno.hexyrs<-pheno.quant %>% group_by(year, cell) %>% 
  summarize(n1=1) %>% group_by(cell) %>% 
  summarize(nyear=sum(n1)) #tally()

#import abundance data
abund<-read_csv(abundance.file) %>%
  filter(cell %in% pheno.hexyrs$cell, ObsMonth %in% c(6:8))


# Figure 1 ----------------------------------------------------------------

####################### FIG 1 A


#summary onset phenology statistics for butterflies overwintering as larvae (RL/BOL)
pheno.ml<-pheno.quant %>%
  filter(code=='RL') %>%
  mutate(ci.weight= (1/(q5_high-q5_low+1))) %>%
  group_by(cell) %>%
  summarize(wmonset=weighted.mean(q5,ci.weight),meanonset=mean(q5, na.rm=T), sdonset=sd(q5, na.rm=T), n=n())


pheno.map.means <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(pheno.ml) # %>%   mutate(n=ifelse(n>17,17,n))

cell.cent<-st_centroid(hex_sf) %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(pheno.ml) #%>%   mutate(n=ifelse(n>17,17,n))

pheno.map.data <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(pheno.hexyrs)

##Create panel A for Figure 1: Data density & mean onset phenology for BOL
fig1a<-ggplot(data = pheno.map.means) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(data=pheno.map.data, color="mediumpurple4", fill=NA) + 
  geom_sf(data=cell.cent, aes(color=wmonset, size=n)) + 
  #geom_sf(aes(fill=wmonset),alpha=0.5) +
  xlim(f1lims.x) + ylim(f1lims.y) + 
  labs(tag="A") + 
  scale_size(range=c(1,8), limits=c(1,18),name="#  Years") + 
  scale_color_viridis(name="Onset DOY",direction=-1) + 
  theme_classic()
fig1a

####################### FIG 1 B

#summary abundance statistics for butterflies overwintering as larvae (RL/BOL)
abund.density<-abund %>%
  filter(group=="RL") %>%
  group_by(cell,Lat, Lng,CountID) %>%
  dplyr::summarize(n=length(unique(ObsYear)), abm=mean(log.abund, na.rm=T))

ab.pts<-st_as_sf(abund.density, coords = c("Lng", "Lat"), crs = 4326) %>%
  st_transform(3857)

fig1b<-ggplot(data = nam_sf) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(data=pheno.map.data, color="mediumpurple4", fill=NA) + 
  geom_sf(data=ab.pts, aes(color=abm, size=n), alpha=0.3) + 
  geom_sf(data=ab.pts, aes(color=abm, size=n),shape=21,  fill=NA) + 
  xlim(f1lims.x) + ylim(f1lims.y) + 
  labs(tag="B") + 
  scale_color_viridis(name="Log abundance", option="A", begin=0.2,end=1,limits = c(min(ab.pts$abm),max(ab.pts$abm))) +
  scale_size(range=c(1,5), limits=c(1,18), name="# Years") + 
  theme_classic()
fig1b

####################### FIG 1 C - phenology across years
load(pheno.dev.file)

pheno.dmean<-pheno.input.dev %>% group_by(code, year) %>% 
  summarize(onset.dev=onset.dev, wm=weighted.mean(onset.dev, onset.ci)) %>% 
  mutate(group=ifelse(code=="RE","BOE",ifelse(code=="RL","BOL","BOP")), year=year+2000) %>%
  filter(year %in% temporal.domain)

fig1c<-ggplot(data=filter(pheno.dmean, year<2018), aes(x=year, y=onset.dev, color=group)) + 
  geom_smooth(method="lm", aes(fill=group),color="black", alpha=0.5,linetype=2) + 
  geom_line(data=pheno.dmean, aes(x=year, y=wm, color=group), size=.8) + 
  scale_color_manual(values= c("BOE" = ows.colors[1],"BOL" = ows.colors[2],"BOP" = ows.colors[3])) + 
  scale_fill_manual(values= c("BOE" = ows.colors[1],"BOL" = ows.colors[2],"BOP" = ows.colors[3])) + 
  xlim(2000,2017) + 
  labs(x="Year", y="Standardized adult onset", color="Group", fill="Group", tag="C") 
fig1c

####################### FIG 1 D - abundance across years
abund.mean<-abund %>% group_by(group, ObsYear) %>% dplyr::summarize(Abundance=log(mean(abund.bph, na.rm=T)), Ab1=exp(Abundance)) %>%
  dplyr::mutate(group=ifelse(group=="RE","BOE",ifelse(group=="RL","BOL","BOP")))


fig1d<-ggplot(data=abund.mean, aes(x=ObsYear, y=Abundance, color=group)) + 
  geom_smooth(method="lm", aes(fill=group),color="black", linetype=2) + 
  geom_line(aes(color=group), size=0.8) + 
  scale_color_manual(values= c("BOE" = ows.colors[1],"BOL" = ows.colors[2],"BOP" = ows.colors[3])) + 
  scale_fill_manual(values= c("BOE" = ows.colors[1],"BOL" = ows.colors[2],"BOP" = ows.colors[3])) + 
  labs(x="Year", y="Log abundance", color="Group", fill="Group", tag="D")
fig1d

#(fig1<-grid.arrange(pheno.meanL2,plot.abund.sites, nrow=1))
(fig1<-grid.arrange(arrangeGrob(fig1a, fig1b, fig1c, fig1d, nrow=2, heights=c(2,1.5))))


ggsave(filename="output/figures/Fig1.202201.png",fig1,width = 12, height = 10)



# Supplemental Figures ----------------------------------------------------



#data density (# years with metrics) by cell & species group (overwinter stage)
pheno.dd.all<-pheno.quant %>%
  group_by(cell, code) %>%
  tally()
  
#summary onset phenology by cell & species group (overwinter stage)
pheno.on.all<-pheno.quant %>%
  dplyr::mutate(ci.weight= (1/(q5_high-q5_low+1))) %>%
  group_by(cell,code) %>%
  dplyr::summarize(wmonset=weighted.mean(q5,ci.weight),meanonset=mean(q5, na.rm=T), sdonset=sd(q5, na.rm=T), n=n())

#xcentroids <- hex_sf %>% st_centroid()  

pheno.map.data2 <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(pheno.dd.all)


pheno.dd<-ggplot(data = pheno.map.data) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=nyear), alpha=0.4)+geom_sf(data=filter(pheno.map.data, cell==766),fill='black') +
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Phenology metric temporal scope") + 
  scale_fill_viridis(name="# Years") + theme_classic()
pheno.dd
ggsave(filename="output/figures/pheno.datamap.png",pheno.dd)

pheno.dd+geom_sf(data=filter(pheno.map.data, cell==766),fill='black')


pheno.meanL<-ggplot(data = pheno.map.means) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=wmonset),alpha=0.5) +
  xlim(f1lims.x) + ylim(f1lims.y) + 
  labs(title="Onset (BOL weighted mean)") + 
  scale_fill_viridis(name="Onset DOY") + 
  scale_alpha(range=rev(c(0.2,0.8)), guide="none") + 
  theme_classic()
pheno.meanL
ggsave(filename="output/figures/pheno.meanO.png",pheno.meanL)

pheno.map.2 <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(pheno.on.all) 


pheno.meanOn<-ggplot(data = pheno.map.2) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=wmonset),alpha=0.5) +
  xlim(f1lims.x) + ylim(f1lims.y) + 
  labs(title="Onset (weighted mean)") + 
  scale_fill_viridis(name="Onset DOY") + 
  scale_alpha(range=rev(c(0.2,0.8)), guide="none") + 
  theme_classic() + facet_wrap(~code, nrow=1)
pheno.meanOn



#### Dev model
load("data/derived/newphenodev.RData")
pheno.map2<-pheno.input.dev %>%
  group_by(cell, code) %>%
  tally(name="nyears")

## Cells data density per code
pheno.map.devs <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(pheno.map2)

pheno.dev.map<-ggplot(data = pheno.map.devs) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=nyears),alpha=0.5) +
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Data density for phenology deviations") + 
  scale_fill_viridis(name="# Years") + 
  scale_alpha(range=rev(c(0.2,0.8)), guide="none") + 
  theme_classic() + facet_wrap(~code)

pheno.grp.map<-ggplot(data = pheno.map.data2) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=n),alpha=0.5) +
  geom_sf(data=pheno.map.devs, color="black", size=1.6, fill=NA) + 
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Data density for annual deviation phenology metrics") + 
  scale_fill_viridis(name="# Years") + 
  scale_alpha(range=rev(c(0.2,0.8)), guide="none") + 
  theme_classic() + facet_wrap(~code)

pheno.grp.map
pheno.dev.map

ggsave(filename="output/figures/Fig.pheno.dd.png",pheno.grp.map,width = 10, height = 5)



### Abundance figures
abund<-read_csv("data/derived/naba_OWS_abundances.csv")
#%>%
#  dplyr::select(cell, year=ObsYear, doy, CountID, SurveyID, ObsMonth, ows.grp=group, abund.bph, log.abund, SR)

abund.density<-abund %>%
  filter(group=="RL") %>%
  group_by(cell,Lat, Lng,CountID) %>%
  dplyr::summarize(n=length(unique(ObsYear)),srm=mean(SR,na.rm=T), abm=mean(log.abund, na.rm=T))

#load("data/abund.input.RData")

ab.pts<-st_as_sf(abund.density, coords = c("Lng", "Lat"), crs = 4326) %>%
  st_transform(3857)


#### FIGURE 1 B

fig1b<-ggplot(data = nam_sf) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(data=pheno.map.data, color="mediumpurple4", fill=NA) + 
  geom_sf(data=ab.pts, aes(color=abm, size=n), alpha=0.3) + 
  geom_sf(data=ab.pts, aes(color=abm, size=n),shape=21,  fill=NA) + 
  xlim(f1lims.x) + ylim(f1lims.y) + 
  labs(tag="B") + 
  scale_color_viridis(name="Log abundance", option="A", begin=0.2,end=1,limits = c(min(ab.pts$abm),max(ab.pts$abm))) +
  #scale_color_viridis(name="Log abundance") + 
  scale_size(range=c(1,5), limits=c(1,18), name="# Years") + 
  #scale_alpha(range=rev(c(0.2,0.8)), guide=F) + 
  theme_classic()
fig1b
#ggsave(filename="output/figures/plot.abund.sites.png",fig1b,width = 7, height = 6,)



############################
## Figure 1C
## 
#load("data/derived/pheno.RData")
pheno.dmean<-pheno.input.dev %>% group_by(code, year) %>% dplyr::summarize(wm=weighted.mean(onset.dev, onset.ci)) %>% dplyr::mutate(group=ifelse(code=="RE","BOE",ifelse(code=="RL","BOL","BOP")), year=year+2000)
pheno.dev<-pheno.input.dev %>% dplyr::mutate(group=ifelse(code=="RE","BOE",ifelse(code=="RL","BOL","BOP")), year=year+2000)

mycolors<-viridis_pal()(6)[c(1,3,6)] #c('#440154','#39568C','#FDE725')  #viridis(20, option="D")[c(1,6,20)]
#mycolors<-c("purple","blue","mustardyellow")
fig1c<-ggplot(data=filter(pheno.dev, year<2018), aes(x=year, y=onset.dev, color=group)) + 
  geom_smooth(method="lm", aes(fill=group),color="black", alpha=0.5,linetype=2) + 
  geom_line(data=pheno.dmean, aes(x=year, y=wm, color=group)) + 
  #scale_color_viridis(discrete=T, begin = 0, end = 0.8, option="magma") + 
  #scale_fill_viridis(discrete=T, begin = 0, end = 0.8, option="magma") + 
  scale_color_manual(values= c("BOE" = mycolors[1],"BOL" = mycolors[2],"BOP" = mycolors[3])) + 
  scale_fill_manual(values= c("BOE" = mycolors[1],"BOL" = mycolors[2],"BOP" = mycolors[3])) + 
  xlim(2000,2017) + 
  labs(x="Year", y="Standardized adult onset", color="Group", fill="Group", tag="C") 
  
fig1c

abund.mean<-abund %>% group_by(group, ObsYear) %>% dplyr::summarize(Abundance=mean(log.abund, na.rm=T), Ab1=exp(Abundance)) %>%
  dplyr::mutate(group=ifelse(group=="RE","BOE",ifelse(group=="RL","BOL","BOP")))


fig1d<-ggplot(data=abund.mean, aes(x=ObsYear, y=Abundance, color=group)) + 
  geom_smooth(method="lm", aes(fill=group),color="black", linetype=2) + 
  geom_line(aes(color=group)) + 
  #scale_color_viridis(discrete=T, begin = 0, end = 0.8, option="plasma") + 
  #scale_fill_viridis(discrete=T, begin = 0, end = 0.8, option="plasma") + 
  scale_color_manual(values= c("BOE" = mycolors[1],"BOL" = mycolors[2],"BOP" = mycolors[3])) + 
  scale_fill_manual(values= c("BOE" = mycolors[1],"BOL" = mycolors[2],"BOP" = mycolors[3])) + 
  labs(x="Year", y="Log abundance", color="Group", fill="Group", tag="D")
fig1d

#(fig1<-grid.arrange(pheno.meanL2,plot.abund.sites, nrow=1))
(fig1<-grid.arrange(arrangeGrob(fig1a, fig1b, fig1c, fig1d, nrow=2, heights=c(2,1.5))))


ggsave(filename="output/figures/Fig1.1213.png",fig1,width = 12, height = 10)



#######################
#Supplemental Figures

green.for<-readRDS("data/envir/MidGreenup-2020-08-06-forest.rds") %>%
  dplyr::select(year, cell, gr_mn_for=gr_mn, gr_pfor=gr_pcell, cell_lat, cell_lng) %>%
  filter(cell %in% STUDYCELLS, !is.na(gr_mn_for)) 

green.for1<-green.for %>%
  group_by(cell) %>%
  summarize(meangreen=mean(gr_mn_for, na.rm=T))

greenup.map.data <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(green.for1)

greenup.map<-ggplot(data = greenup.map.data) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=meangreen),alpha=0.5) +
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Mean greenup DOY") + 
  scale_fill_viridis(name="Forest Greenup DOY", direction=-1) + 
  scale_alpha(range=rev(c(0.2,0.8)), guide="none") + 
  theme_classic()
greenup.map





##Lag map
green.all<-readRDS("data/envir/MidGreenup-2020-08-06-all.rds") %>% 
  dplyr::select(year, cell, gr_mn_all=gr_mn, gr_pall=gr_pcell, cell_lat, cell_lng) %>%
  filter(cell %in% STUDYCELLS, !is.na(gr_mn_all))


#Greenup: calculate half-max DOY of open pixels
greendev<-merge(green.all,green.for, by=intersect(names(green.all),names(green.for)), all.x=T, all.y=T) %>%
  mutate(p_open=gr_pall-gr_pfor) %>%
  mutate(gr_mn_open=ifelse(p_open>0,((gr_mn_all*gr_pall - gr_mn_for*gr_pfor)/(gr_pall-gr_pfor)),NA), gr_mn_lag=(gr_mn_open-gr_mn_for)) %>%
  group_by(cell) %>%
  summarize(gr_mn_open=mean(gr_mn_open, na.rm=T), meanlag=mean(gr_mn_lag, na.rm=T))



greenuplag.map.data <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(greendev)

greenuplag.map<-ggplot(data = greenuplag.map.data) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=meanlag),alpha=0.5) +
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Mean lags to open canopy greenup") + 
  scale_fill_viridis(name="Mean lag") + 
  scale_alpha(range=rev(c(0.2,0.8)), guide="none") + 
  theme_classic()
greenuplag.map


greenup.open.map<-ggplot(data = greenuplag.map.data) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=gr_mn_open),alpha=0.5) +
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Open canopy greenup") + 
  scale_fill_viridis(name="Mean DOY", direction=-1) + 
  scale_alpha(range=rev(c(0.2,0.8)), guide="none") + 
  theme_classic()
greenup.open.map

env2<-read_csv("data/derived/envir2.csv")



env2a<-env2 %>%
  group_by(cell) %>%
  summarize(spring.gdd=mean(spring.gdd, na.rm=T), pc1=mean(pc1, na.rm=T), pc2=mean(pc2, na.rm=T))
env2b<-env2 %>%
  group_by(year) %>%
  summarize(spring.gdd=mean(spring.gdd, na.rm=T), pc1=mean(pc1, na.rm=T), pc2=mean(pc2, na.rm=T))

#### GDD 
gdd.map.data <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(env2a)

gdd.map<-ggplot(data = gdd.map.data) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=spring.gdd),alpha=0.5) +
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Mean Spring GDD") + 
  scale_fill_viridis(name="GDD") + 
  scale_alpha(range=rev(c(0.2,0.8)), guide="none") + 
  theme_classic()
gdd.map


# Spatiotemporal PC maps (Supplemental figure 3.3)
stpc.data<-read_csv("data/derived/spatemp_env.csv") %>%
  group_by(cell) %>%
  summarize(ST.PC1=mean(ST.PC1, na.rm=T),ST.PC2=mean(ST.PC2, na.rm=T))

stpc.map.data <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(stpc.data)

stpc1.map<-ggplot(data = stpc.map.data) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=ST.PC1),alpha=0.5) +
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Mean ST.PC1") + 
  scale_fill_viridis(name="ST.PC1") + 
  theme_classic()
stpc1.map

stpc2.map<-ggplot(data = stpc.map.data) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=ST.PC2),alpha=0.5) +
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Mean ST.PC2") + 
  scale_fill_viridis(name="ST.PC2") + 
  theme_classic()
stpc2.map


(figs33<-grid.arrange(arrangeGrob(stpc1.map, stpc2.map, nrow=1)))


ggsave(filename="output/figures/FigS3.3.png",figs33,width = 10, height = 5)




