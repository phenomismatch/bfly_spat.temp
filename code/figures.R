###FIGURES
library(tidyverse)
library(sf)

theme_set(theme_classic(base_size = 15))


#phenology data
pheno.datafile<-"data/derived/adult_bfly_phenometrics_noCountCircles.csv"
pheno.orig<-read_csv(pheno.datafile) %>% rename(cell=HEXcell)
hexyrs<-pheno.orig %>% dplyr::select(year, cell)  %>% group_by(year, cell) %>% dplyr::summarize(n1=1) %>% group_by(cell) %>% summarize(nyear=sum(n1)) #tally()
nrow(pheno.orig)
#pheno.devs<-pheno.input.dev %>%
#  group_by(cell, code) %>%
#  summarize(meanonset=mean(onset, na.rm=T), sdonset=sd(onset, na.rm=T), nyears=n())

pheno.ml<-pheno.orig %>%
  filter(code=='RL') %>%
  mutate(ci.weight= (1/(q5_high-q5_low+1))) %>%
  group_by(cell) %>%
  summarize(wmonset=weighted.mean(q5,ci.weight),meanonset=mean(q5, na.rm=T), sdonset=sd(q5, na.rm=T))


hex_sf.ll <- read_sf("data/maps/hex_grid_crop.shp") 

hex_sf<-st_transform(hex_sf.ll, 3857) #3857 is projected pseudo-mercator, 4326 is lat/long


# check if centroid is in polygon
xcentroids <- hex_sf %>% st_centroid()  
llcentroids<-st_transform(xcentroids, 4326) %>%
  st_crop(xmin = -95, ymin = 35, xmax = -60, ymax = 50)


STUDYCELLS<-as.numeric(llcentroids$cell)
save(STUDYCELLS,file="data/spatial.domain.RData")

hex_sf<-hex_sf %>% filter(cell %in% STUDYCELLS)
summary(hex_sf)


nam_sf <- read_sf("data/maps/ne_50m_admin_1_states_provinces_lakes.shp") %>%
  st_crop(xmin = -95, ymin = 25, xmax = -60, ymax = 50)

nam_sf<-st_transform(nam_sf, 3857)
nam_sf<-nam_sf %>%
  filter(sr_adm0_a3 %in% c("CAN", "USA"), iso_3166_2 != "US-AK", iso_3166_2 != "US-HI") 

library("rnaturalearth")
library("rnaturalearthdata")
library("rgeos")

spdf.n.am <- ne_countries(country = c("united states of america","canada","mexico"), returnclass = 'sf') 
basmap<-st_transform(spdf.n.am,  "+init=epsg:3857")
pt1 = st_point(c(-100,24))
pt2 = st_point(c(-60,50))

bounds<-st_sfc(pt1, pt2) # %>% st_transform(crs ="+init=epsg:3857")
sf_points <- 
  bounds %>% 
  st_as_sf(coords = c('a', 'b'), crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs ="+init=epsg:3857")
box1<-st_bbox(sf_points)




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

pheno.meanL<-ggplot(data = pheno.map.means) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=wmonset),alpha=0.5) +
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Onset (weighted mean)", subtitle="Sp. overwintering as caterpillars") + 
  scale_fill_viridis(name="Onset DOY") + 
  scale_alpha(range=rev(c(0.2,0.8)), guide="none") + 
  theme_classic()
pheno.meanL
ggsave(filename="output/figures/pheno.meanO.png",pheno.meanL)

pheno.devs<-ggplot(data = pheno.map.devs) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=nyears),alpha=0.5) +
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Data density for temporal deviations") + 
  scale_fill_viridis(name="# Years") + 
  scale_alpha(range=rev(c(0.2,0.8)), guide="none") + 
  theme_classic() + facet_wrap(~code)
pheno.devs
ggsave(filename="output/figures/pheno.devs.png",pheno.devs)

pheno.map1<-pheno.quant %>%
  group_by(cell, code) %>%
  tally(name="nyears")

## Cells data density per code
pheno.map.data <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(pheno.map1)

pheno.devs<-ggplot(data = pheno.map.data) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(aes(fill=nyears),alpha=0.5) +
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Data density for phenometrics") + 
  scale_fill_viridis(name="# Years") + 
  scale_alpha(range=rev(c(0.2,0.8)), guide="none") + 
  theme_classic() + facet_wrap(~code)

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
  filter(group=="RL", cell %in% hexyrs$cell) %>%
  group_by(Lat, Lng) %>%
  summarize(n=n(), srm=mean(SR,na.rm=T), abm=mean(log.abund, na.rm=T))

load("data/abund.input.RData")

ab.pts<-st_as_sf(abund.density, coords = c("Lng", "Lat"), crs = 4326) %>%
  st_transform(3857)


plot.abund.sites<-ggplot(data = nam_sf) +
  geom_sf(data=nam_sf, color="darkgray", fill="white") + 
  geom_sf(data=pheno.map.data, color="mediumpurple4", fill=NA) + 
  geom_sf(data=ab.pts, aes(color=abm, size=n)) + 
  xlim(box1[c(1,3)]) + ylim(box1[c(2,4)]) + 
  labs(title="Abundance surveys", subtitle="Mean abundance for sp. overwintering as larvae") + 
  scale_color_viridis(name="Log abundance") + 
  scale_size(range=c(1,4), name="# Survey Years") + 
  #scale_alpha(range=rev(c(0.2,0.8)), guide=F) + 
  theme_classic()
plot.abund.sites
ggsave(filename="output/figures/plot.abund.sites.png",plot.abund.sites,width = 7, height = 6,)

abund.best<-lmer(log.abund~-1+ows.grp+logab.py+abslag+ows.grp:warmearly+ows.grp:warmlateopen+ows.grp:on.dev+ows.grp:year+ows.grp:as.factor(ObsMonth)+ows.grp:doy+FR.dev+(1|cell) + (1|CountID:cell), data=naba.1)
c1<- predict(abund.best, re.form=~0)
#ab.pred <- cbind(naba.1, predict(abund.best)) #, interval = 'confidence'))

ggplot(data=ab.pred, aes(x=on.dev, y=year, color=logabund))


library(ggeffects)
library(lme4)

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
