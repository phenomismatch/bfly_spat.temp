###FIGURES
library(ggeffects)
library(lme4)
library(sf)
library(viridis)
library(tidyverse)
library(gridExtra)
theme_set(theme_classic(base_size = 15))

load("data/spatial.domain.RData")
#phenology data
pheno.datafile<-"data/derived/adult_bfly_phenometrics_noCountCircles_withFull2020Data.csv"
pheno.quant<-read_csv(pheno.datafile) %>% dplyr::rename(cell=HEXcell) %>%
  filter(cell %in% STUDYCELLS, between(q50,152,243), !is.na(q5))

hexyrs<-pheno.quant %>% dplyr::select(year, cell)  %>% group_by(year, cell) %>% dplyr::summarize(n1=1) %>% group_by(cell) %>% dplyr::summarize(nyear=sum(n1)) #tally()
#nrow(pheno.quant)
#pheno.devs<-pheno.input.dev %>%
#  group_by(cell, code) %>%
#  summarize(meanonset=mean(onset, na.rm=T), sdonset=sd(onset, na.rm=T), nyears=n())

pheno.ml<-pheno.quant %>%
  filter(code=='RL') %>%
  dplyr::mutate(ci.weight= (1/(q5_high-q5_low+1))) %>%
  group_by(cell) %>%
  dplyr::summarize(wmonset=weighted.mean(q5,ci.weight),meanonset=mean(q5, na.rm=T), sdonset=sd(q5, na.rm=T), n=n())


hex_sf.ll <- read_sf("data/maps/hex_grid_crop.shp") 

hex_sf<-st_transform(hex_sf.ll, 3857) #3857 is projected pseudo-mercator, 4326 is lat/long


# check if centroid is in polygon
xcentroids <- hex_sf %>% st_centroid()  
#llcentroids<-st_transform(xcentroids, 4326) %>%
#  st_crop(xmin = -95, ymin = 35, xmax = -60, ymax = 50)


#STUDYCELLS<-as.numeric(llcentroids$cell)
#save(STUDYCELLS,file="data/spatial.domain.RData")

hex_sf<-hex_sf %>% filter(cell %in% STUDYCELLS)
summary(hex_sf)


nam_sf <- read_sf("data/maps/ne_50m_admin_1_states_provinces_lakes.shp") %>%
  st_crop(xmin = -95, ymin = 25, xmax = -60, ymax = 50)

nam_sf<-st_transform(nam_sf, 3857)
nam_sf<-nam_sf %>%
  filter(sr_adm0_a3 %in% c("CAN", "USA"), iso_3166_2 != "US-AK", iso_3166_2 != "US-HI") 

#library("rnaturalearth")
#library("rnaturalearthdata")
#library("rgeos")

#spdf.n.am <- ne_countries(country = c("united states of america","canada","mexico"), returnclass = 'sf') 
#basmap<-st_transform(spdf.n.am,  "+init=epsg:3857")

pt1 = st_point(c(-100,24))
pt2 = st_point(c(-60,50))

bounds<-st_sfc(pt1, pt2) # %>% st_transform(crs ="+init=epsg:3857")
sf_points <- 
  bounds %>% 
  st_as_sf(coords = c('a', 'b'), crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs ="+init=epsg:3857")
box1<-st_bbox(sf_points)


cell.cent<-xcentroids %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(pheno.ml)

pheno.map.data <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(hexyrs)

pheno.map.means <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(pheno.ml)

#pheno.map.devs <- hex_sf %>%
#  mutate_at(c("cell"), ~as.numeric(.)) %>%
#  right_join(pheno.dev)


pheno.dd<-ggplot(data = pheno.map.data) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=nyear), alpha=0.4) +
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Phenology metric temporal scope") + 
  scale_fill_viridis(name="# Years") + theme_classic()
pheno.dd
ggsave(filename="output/figures/pheno.datamap.png",pheno.dd)

f1lims.x<-c(-10700000,-6679169)
f1lims.y<-c(3500000,6500000)

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


####################### FIG 1 A
fig1a<-ggplot(data = pheno.map.means) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(data=pheno.map.data, color="mediumpurple4", fill=NA) + 
  geom_sf(data=cell.cent, aes(color=wmonset, size=n)) + 
  #geom_sf(aes(fill=wmonset),alpha=0.5) +
  xlim(f1lims.x) + ylim(f1lims.y) + 
  labs(tag="A") + 
  scale_size(range=c(3,7), limits=c(1,18),name="#  Years") + 
  scale_color_viridis(name="Onset DOY",direction=-1) + 
  theme_classic()
fig1a



#### Dev model

pheno.map2<-pheno.dev %>%
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

pheno.dev.map



### Abundance figures
abund<-read_csv("data/derived/naba_OWS_abundances.csv")
#%>%
#  dplyr::select(cell, year=ObsYear, doy, CountID, SurveyID, ObsMonth, ows.grp=group, abund.bph, log.abund, SR)

abund.density<-abund %>%
  filter(group=="RL", cell %in% hexyrs$cell, ObsMonth %in% c(6:8)) %>%
  group_by(cell,Lat, Lng,CountID) %>%
  dplyr::summarize(n=length(unique(ObsYear)),srm=mean(SR,na.rm=T), abm=mean(log.abund, na.rm=T))

#load("data/abund.input.RData")

ab.pts<-st_as_sf(abund.density, coords = c("Lng", "Lat"), crs = 4326) %>%
  st_transform(3857)


#### FIGURE 1 B

fig1b<-ggplot(data = nam_sf) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(data=pheno.map.data, color="mediumpurple4", fill=NA) + 
  geom_sf(data=ab.pts, aes(color=abm, size=n), alpha=0.15) + 
  geom_sf(data=ab.pts, aes(color=abm, size=n),shape=21,  fill=NA) + 
  xlim(f1lims.x) + ylim(f1lims.y) + 
  labs(tag="B") + 
  scale_color_viridis(name="Log abundance", option="A", begin=0.2,end=1,limits = c(min(ab.pts$abm),max(ab.pts$abm))) +
  #scale_color_viridis(name="Log abundance") + 
  scale_size(range=c(1,4), limits=c(1,18), name="# Years") + 
  #scale_alpha(range=rev(c(0.2,0.8)), guide=F) + 
  theme_classic()
fig1b
#ggsave(filename="output/figures/plot.abund.sites.png",fig1b,width = 7, height = 6,)



############################
## Figure 1C
## 
load("data/derived/pheno.RData")
pheno.dmean<-pheno.dev %>% group_by(code, year) %>% dplyr::summarize(wm=weighted.mean(onset.dev, onset.ci)) %>% dplyr::mutate(group=ifelse(code=="RE","BOE",ifelse(code=="RL","BOL","BOP")))
pheno.dev<-pheno.dev %>% dplyr::mutate(group=ifelse(code=="RE","BOE",ifelse(code=="RL","BOL","BOP")))

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

abund.mean<-abund %>% group_by(group, ObsYear) %>% dplyr::summarize(Abundance=mean(log.abund, na.rm=T)) %>%
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


ggsave(filename="output/figures/Fig1.1027.png",fig1,width = 12, height = 10)





abund.best<-lmer(log.abund~-1+ows.grp+logab.py+abslag+ows.grp:warmearly+ows.grp:warmlateopen+ows.grp:on.dev+ows.grp:year+ows.grp:as.factor(ObsMonth)+ows.grp:doy+FR.dev+(1|cell) + (1|CountID:cell), data=naba.1)
c1<- predict(abund.best, re.form=~0)
#ab.pred <- cbind(naba.1, predict(abund.best)) #, interval = 'confidence'))

ggplot(data=ab.pred, aes(x=on.dev, y=year, color=logabund))



pr <- ggpredict(abund.best, c("on.dev","ows.grp"))
pr
plot(pr) + scale_color_viridis_d()

pr <- ggpredict(abund.best, c("warmearly","ows.grp"))
pr
plot(pr) + scale_color_viridis_d()


pr <- ggpredict(abund.best, c("doy","ows.grp","ObsMonth"))
pr
plot(pr) + scale_color_viridis_d()


##Plot environment drivers of abundance 
env.var<-list("FR.dev","FFD.dev","warmearly","warmlateopen")
grefx<-function(e.var) {
  grep(e.var,row.names(summary(abund.best)$coefficients))
}
allvar<-unlist(lapply(env.var, grefx))

model.env.coefs<-as.data.frame(summary(abund.best)$coefficients[allvar,]) 
mec<-cbind(model.env.coefs,as.data.frame(matrix(unlist(str_split(row.names(model.env.coefs),":")),byrow=T,ncol=2)))
names(mec)<-c(names(model.env.coefs),"ows","env.var")                
mec<-mec %>% dplyr::mutate(ows=ifelse(ows=="ows.grpRE","RE",ifelse(ows=="ows.grpRL","RL","RP")),
                           axisval=as.numeric(factor(env.var))+ifelse(ows=="RE",0.2,ifelse(ows=="RL",0.1,0)),
                           lowci=Estimate-2*`Std. Error`,highci=Estimate+2*`Std. Error`,
                           xval=as.numeric(factor(env.var)),
                           sig=ifelse(`Pr(>|t|)`<0.05,"*",""))
mec2<-mec %>% group_by(env.var,xval) %>% tally()
ggplot(data=mec, aes(x=axisval,y=Estimate, ymin = lowci, ymax = highci, color=ows)) + 
  geom_pointrange() + coord_flip() + 
  scale_x_continuous(breaks=c(mec2$xval),labels=c("Cold days deviation","PC1 (warm, early)","PC2 (warm; late open greenup)"),name="") + 
  scale_color_viridis_d() + theme_classic() + geom_hline(yintercept=0, linetype="dashed") + 
  labs(title="Environmental factors in abundance model", y="Coefficient estimates") + 
  annotate("text", x = (mec$axisval+0.1), y = mec$Estimate, label = mec$sig)
  


## Results tables

model.stats<- function(finalmodel) {
  model.coefs<-as.data.frame(summary(finalmodel)$coefficients) %>%
    dplyr::mutate(Estimate=round(Estimate,3),sig=ifelse(`Pr(>|t|)`<0.05,1,0)) %>%
    dplyr::select(Estimate, `Pr(>|t|)`,sig)
  model.r2<-as.data.frame(cbind(round(c(r.squaredGLMM(finalmodel)),3),c(rep(NA,2)),c(rep(NA,2))))
  row.names(model.r2)<-c("r2m","r2c")
  names(model.r2)<-names(model.coefs)
  model.table<-rbind(model.coefs, model.r2)
  
  }



onset.dev.best<-lmer(dev_onset~-1+code+code:year+warmearly+(1|cell), data=pheno.input.dev)
r.squaredGLMM(onset.dev.best)     

onset.best<-lmer(onset~-1+code+code:year+warmearly+(1|cell), data=pheno.input)
r.squaredGLMM(onset.best)     


abund.result<-model.stats(abund.best)
onset.result<-model.stats(onset.best)
onset.dev.result<-model.stats(onset.dev.best)
write.csv(abund.result, file="output/abund.model.csv")
write.csv(onset.result, file="output/onset.model.csv")
write.csv(onset.dev.result, file="output/onset.dev.model.csv")



### use magma or plasma for fig 1c, 1d