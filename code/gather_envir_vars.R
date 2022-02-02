#Environmental variables for phenology & abundance modeling
#E Larsen, Georgetown U, Updated 2022-01

##Sources
#GDD, frost-free days and frozen days calcs use CPC data:
#Cite: CPC Global Temperature data provided by the NOAA/OAR/ESRL PSL, Boulder, Colorado, USA, from their Web site at https://psl.noaa.gov/data/gridded/data.cpc.globaltemp.html

#libraries
library(tidyverse)
library(ggcorrplot)
library(gridExtra)
library(sp)
library(viridis)
#Set Parameters
rundat<-Sys.Date()
#input files
clim.csv<-"data/derived/therm_metrics.csv"
forest.greenup.input<-"data/envir/MidGreenup-2021-11-30-forest.rds"
open.greenup.input<-"data/envir/MidGreenup-2021-11-30-open.rds"
#output files
env.dev.csv<-"data/derived/envDevs.csv"
env.st.csv<-"data/derived/spatemp_env.csv"
eq_doys<-c(80,172,266,356) #DOYs for equinox & solstice dates
#cells in study area
load("data/spatial.domain.RData")


## SUMMER CLIMATE
#Calculate summary seasonal GDD values for each hex cell & year:
## Spring GDD (from DOY 80 to DOY 172: approximating spring equinox to summer solstice)

#CPC data has spring gdd, summer gdd, winter cold days (tmax<=0), and winter warm days (tmin>0)
clim.data<-read_csv(clim.csv) %>%
  filter(cell %in% STUDYCELLS)

## GREENUP
###Timing of half-max greenup (forest pixels)

#Greenup: halfmax DOY forest pixels
green.for<-readRDS(forest.greenup.input) %>%
  dplyr::select(year, cell, gr_mn_for=gr_mn, gr_pfor=gr_pcell, cell_lat, cell_lng) %>%
  filter(cell %in% STUDYCELLS, !is.na(gr_mn_for))

#Greenup: half-max DOY all pixels
green.open<-readRDS(open.greenup.input) %>% 
  dplyr::select(year, cell, gr_mn_open=gr_mn, gr_popen=gr_pcell, cell_lat, cell_lng) %>%
  filter(cell %in% STUDYCELLS, !is.na(gr_mn_open))

#Greenup: calculate greenup lag between forest and open pixels
greendev<-merge(green.open,green.for, by=intersect(names(green.open),names(green.for)), all.x=T, all.y=T) %>%
  mutate(gr_mn_lag=(gr_mn_open-gr_mn_for))

#visualize greenup data
ggplot(data=filter(greendev, year>2009), aes(x=gr_mn_for, y=gr_mn_lag, color=cell_lat)) + geom_point() + 
  labs(x="Greenup in forest pixels", y="Lag of open pixels", color="latitude") + 
  theme_minimal()
#Combined dataframe for environmental variables
env_var<-merge(greendev, clim.data, by=intersect(names(greendev), names(clim.data)), all.x=T, all.y=T)

# Deviation metrics -------------------------------------------------------

env_baseline<-env_var %>%
  mutate(base.year=ifelse(year>2015,1,NA)) %>%
  mutate(g.lag=gr_mn_lag*base.year, g.for=gr_mn_for*base.year,
         cold=colddays*base.year, warm=warmdays*base.year, 
         spring=spring.gdd*base.year,summer=summer.gdd*base.year) %>%
  group_by(cell) %>%
  mutate(base.for=mean(g.for, na.rm=T), base.lag=mean(g.lag, na.rm=T),
         base.cold=mean(cold, na.rm=T), base.warm=mean(warm, na.rm=T), 
         base.spring=mean(spring, na.rm=T), base.summer=mean(summer, na.rm=T)) 
  
env_deviations<-env_var %>%
  filter(year<2020) %>%
  mutate(base.year=ifelse(year>2015,1,NA)) %>%
  mutate(g.lag=gr_mn_lag*base.year, g.for=gr_mn_for*base.year,
         cold=colddays*base.year, warm=warmdays*base.year, 
         spring=spring.gdd*base.year,summer=summer.gdd*base.year) %>%
  group_by(cell) %>%
  mutate(base.for=mean(g.for, na.rm=T), base.lag=mean(g.lag, na.rm=T),
            base.cold=mean(cold, na.rm=T), base.warm=mean(warm, na.rm=T), 
            base.spring=mean(spring, na.rm=T), base.summer=mean(summer, na.rm=T)) %>%
  mutate(for.green.dev=gr_mn_for-base.for, open.lag.dev=gr_mn_lag-base.lag,
         cold.dev=colddays-base.cold, warm.dev=warmdays-base.warm,
         spring.gdd.dev=spring.gdd-base.spring, summer.gdd.dev=summer.gdd-base.summer) %>%
  select(year:cell_lng,for.green.dev:summer.gdd.dev)


##PCA for correlated variables
dev.corr <- round(cor(na.omit(env_deviations %>% select(for.green.dev:summer.gdd.dev))), 2)
ggcorrplot(dev.corr, type="lower",lab = TRUE)
#forest greenup deviation and spring gdd deviation are negatively correlated. 
#Winter warm days & cold days deviations are negatively correlated, but never used in the same model.

env.dev.pc<-prcomp(na.omit(env_deviations[,c("spring.gdd.dev","for.green.dev")]), scale=T)
#scale the principal components such that the pcs are scaled to days of greenup shift
x1<-as.data.frame(env.dev.pc$x)*-1*(env.dev.pc$scale[2]/abs(env.dev.pc$rotation[2,1]))
names(x1)<-c("dev.pc1","dev.pc2")
env.dev.pc$rotation

##Visualize PCA
pca.viz<-bind_cols(na.omit(env_deviations[,c("spring.gdd.dev","for.green.dev", "cell_lat", "year")]),x1)

(biplot1<-ggplot(data=na.omit(pca.viz), aes(x=dev.pc1, y=dev.pc2, color=cell_lat)) + geom_point()  + 
    scale_color_viridis() + 
    geom_segment(data=as.data.frame(env.dev.pc$rotation)*env.dev.pc$scale[2]*2, inherit.aes=FALSE, aes(x=0, xend=PC1*2, y=0, yend=PC2*2), color="black", arrow = arrow(angle=30, length=unit(0.2,"cm"),type="open"))  + 
    geom_text(data=as.data.frame(env.dev.pc$rotation)*env.dev.pc$scale[2]*2, inherit.aes=FALSE, aes(x=PC1*2, y=PC2*2, label=c("Greenup deviation","Spring GDD deviation")), size=3.5, hjust=0.5) + 
    labs(x="PC 1", y="PC 2") + theme_classic() )

#Add PCA to environmental deviations dataframe (removes rows with NA greenup or spring gdd)
envir.dev<-bind_cols(env_deviations[which(!is.na(env_deviations$for.green.dev) & !is.na(env_deviations$spring.gdd.dev)),],x1)

##### SAVE ENVIRONMENTAL VARIABLES (ANNUAL DEVIATIONS) FOR MAIN ANALYSIS
write.csv(envir.dev, file=env.dev.csv)

#Some visualizations
ggplot(data=envir.dev, aes(x=spring.gdd.dev, y=for.green.dev, color=cell_lat)) + geom_point()
ggplot(data=envir.dev, aes(x=dev.pc1, y=dev.pc2, color=spring.gdd.dev)) + geom_point() + labs(x="PC1", y="PC2") + theme_minimal()
ggplot(data=envir.dev, aes(x=dev.pc1, y=dev.pc2, color=for.green.dev)) + geom_point() + labs(x="PC1", y="PC2") + theme_minimal()
ggplot(data=envir.dev, aes(x=dev.pc1, y=dev.pc2, color=cell_lat)) + geom_point() + labs(x="PC1", y="PC2") + theme_minimal()

rm(x1)

# Spatiotemporal environmental variables ----------------------------------
summary(env_var)

(figs31a<-ggplot(data=env_var, aes(x=spring.gdd, y=gr_mn_for, color=cell_lat)) + 
  geom_point() + theme_classic() + labs(x="Spring GDD", y="Forest greenup DOY", color="Latitude") +
  scale_color_viridis() + theme(legend.position="right"))

spatemp.cor<-(cor(na.omit(env_var) %>% select(spring.gdd, forest.greenup=gr_mn_for, open.greenup.lag=gr_mn_lag )))
(figs31b<-corrplot::corrplot(spatemp.cor, method="number", type="upper"))
                  
library(gridGraphics)
library(grid)
grid.echo()
P1 <- grid.grab()

(figS3.1<-grid.arrange(arrangeGrob( P1,figs31a, nrow=1)))


ggsave(filename="output/figures/FigS31.png",figS3.1,width = 8, height = 4)

env_var<-env_var %>% 
  rename(forest.greenup=gr_mn_for)
#Principal components of spatiotemporal spring.gdd, forest.greenup
env.spatemp.pc<-prcomp(na.omit(env_var[,c("spring.gdd","forest.greenup")]), scale=T)
#scale the principal components such that the pcs are scaled to days of greenup shift
x1<-as.data.frame(env.spatemp.pc$x)*(env.spatemp.pc$scale[2]/abs((env.spatemp.pc$rotation[2,1])))
names(x1)<-c("ST.PC1","ST.PC2")
#if(sign(env.spatemp.pc$rotation[1,1]<0)) {env.spatemp.pc$rotation[,1]<-env.spatemp.pc$rotation[,1]*-1}
env.spatemp.pc$rotation


##Visualize PCA
st.pca.viz<-bind_cols(na.omit(env_var[,c("spring.gdd","forest.greenup", "cell_lat", "year")]),x1)

(figs32a<-ggplot(data=na.omit(st.pca.viz), aes(x=ST.PC1, y=ST.PC2, color=cell_lat)) + geom_point()  + 
    scale_color_viridis() + xlim(-60,85) +
    geom_segment(data=as.data.frame(env.spatemp.pc$rotation)*env.spatemp.pc$scale[2]*2, inherit.aes=FALSE, aes(x=0, xend=PC1*2, y=0, yend=PC2*2), color="black", arrow = arrow(angle=30, length=unit(0.2,"cm"),type="open"))  + 
    geom_text(data=as.data.frame(env.spatemp.pc$rotation)*env.spatemp.pc$scale[2]*2, inherit.aes=FALSE, aes(x=PC1*2, y=PC2*2, label=c("Spring GDD","Forest greenup")), size=3.5, hjust=0.5) + 
    labs(x="ST.PC 1", y="ST.PC 2", color="Latitude") + theme_classic() )

st.pca.means<-na.omit(st.pca.viz) %>% group_by(year) %>% 
  summarize(`ST.PC1`=mean(ST.PC1, na.rm=T),`ST.PC2`=mean(ST.PC2, na.rm=T)) %>%
  pivot_longer(cols=c(`ST.PC1`,`ST.PC2`), names_to=c("PC"))

(figs32b<-ggplot(data=filter(st.pca.means, year %in% c(2001:2017)), aes(x=year, y=value)) +  
    geom_line(aes(color=PC)) + 
    scale_color_viridis_d(alpha=1, begin=0, end=.6, option="A", guide="none") +
    geom_text(data=filter(st.pca.means, year==2003), hjust=-0.1, aes(y=value*.3+10,color=PC, label=PC)) + 
    #annotate("text",x=c(2004,2004), y=c(11,8), labels=c("ST.PC1", "ST.PC2"), aes(color=stpca.col)) +
    labs(x="Year", y="Mean Spatiotemporal PC value") + theme_classic() )

(figS3.2<-grid.arrange(arrangeGrob(figs32a,figs32b, nrow=1)))

ggsave(filename="output/figures/FigS32.png",figS3.2,width = 10, height = 4)




#Add PCA to environmental deviations dataframe (removes rows with NA greenup or spring gdd)
envir.spatemp<-bind_cols(env_var[which(!is.na(env_var$forest.greenup) & !is.na(env_var$spring.gdd)),],x1)

write.csv(envir.spatemp, file=env.st.csv)
